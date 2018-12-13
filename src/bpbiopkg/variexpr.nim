import ./duko
import ./pedfile
import ./duko
import strutils
import tables
import hts/vcf
import os
import times
import strformat
import docopt


type ISample = object
  i: int
  duk: Duko

type Trio = array[3, ISample] ## kid, dad, mom

type TrioEvaluator* = ref object
  ctx: DTContext
  trios: seq[Trio]
  INFO: Duko
  expressions: seq[Dukexpr]
  names: seq[string]

proc fill[T: int8 | int32 | float32 | string](sample:ISample, name:string, values:var seq[T], nper:int) {.inline.} =
  if nper > 2: return
  if nper == 1:
    sample.duk[name] = values[sample.i]
  else:
    sample.duk[name] = values[(nper*sample.i)..<(nper*(sample.i+1))]

proc fill[T: int8 | int32 | float32 | string](trio:Trio, name:string, values:var seq[T], nper:int) {.inline.} =
  var filled = newSeq[bool](int(values.len / nper))
  for s in trio:
    if filled[s.i]:
        continue
    filled[s.i] = true
    s.fill(name, values, nper)

proc newEvaluator*(kids: seq[Sample], expression: Table[string, string]): TrioEvaluator =
  ## make a new evaluation context for the given string
  var my_fatal: duk_fatal_function = (proc (udata: pointer, msg:cstring) {.stdcall.} =
    stderr.write_line "varianteval fatal error:"
    quit $msg
  )

  result = TrioEvaluator(ctx:duk_create_heap(nil, nil, nil, nil, my_fatal))
  result.ctx.duk_require_stack_top(50000)
  for k, v in expression:
    result.expressions.add(result.ctx.compile(v))
    result.names.add(k)
  for kid in kids:
      result.trios.add([ISample(i:kid.i, duk:result.ctx.newObject(kid.id)),
                        ISample(i:kid.dad.i, duk:result.ctx.newObject(kid.dad.id)),
                        ISample(i:kid.mom.i, duk:result.ctx.newObject(kid.mom.id))])
  result.INFO = result.ctx.newObject("INFO")

proc clear*(ctx:var TrioEvaluator) {.inline.} =
  for trio in ctx.trios.mitems:
    trio[0].duk.clear()
    trio[1].duk.clear()
    trio[2].duk.clear()
  ctx.INFO.clear()

proc set_format_field(ctx: TrioEvaluator, f:FormatField, fmt:FORMAT, ints: var seq[int32], floats: var seq[float32]) =

  if f.vtype == BCF_TYPE.FLOAT:
    if fmt.get(f.name, floats) != Status.OK:
      quit "couldn't get format field:" & f.name
    for trio in ctx.trios:
      trio.fill(f.name, floats, f.n_per_sample)
  elif f.vtype == BCF_TYPE.CHAR:
    discard
  elif f.vtype in {BCF_TYPE.INT32, BCF_TYPE.INT16, BCF_TYPE.INT8}:
    if fmt.get(f.name, ints) != Status.OK:
      quit "couldn't get format field:" & f.name
    for trio in ctx.trios:
      trio.fill(f.name, ints, f.n_per_sample)
  else:
    quit "Unknown field type:" & $f.vtype & " in field:" & f.name

proc set_infos(ctx:var TrioEvaluator, variant:Variant, ints: var seq[int32], floats: var seq[float32]) =
  var istr: string = ""
  var info = variant.info
  for field in info.fields:
    if field.vtype == BCF_TYPE.FLOAT:
      if info.get(field.name, floats) != Status.OK:
        quit "couldn't get field:" & field.name
      if field.n == 1:
          ctx.INFO[field.name] = floats[0]
      else:
          ctx.INFO[field.name] = floats
    elif field.vtype == BCF_TYPE.CHAR:
      when true:
        var ret = info.get(field.name, istr)
        if ret != Status.OK:
          quit "couldn't get field:" & field.name & " status:" & $ret
        # NOTE: all set as a single string for now.
        ctx.INFO[field.name] = $istr
    elif field.vtype in {BCF_TYPE.INT32, BCF_TYPE.INT16, BCF_TYPE.INT8}:
      if info.get(field.name, ints) != Status.OK:
        quit "couldn't get field:" & field.name
      if field.n == 1:
          ctx.INFO[field.name] = ints[0]
      else:
          ctx.INFO[field.name] = ints


proc evaluate*(ctx:var TrioEvaluator, variant:Variant, samples:seq[string]): TableRef[string, seq[string]] =
  ctx.clear()
  var ints = newSeq[int32](3 * variant.n_samples)
  var floats = newSeq[float32](3 * variant.n_samples)

  ## the most expensive part is pulling out the format fields so we pull all fields
  ## and set values for all samples in the trio list.
  ## once all that is done, we evaluate the expressions.
  result = newTable[string, seq[string]]()

  ctx.set_infos(variant, ints, floats)

  # file the format fields
  var fmt = variant.format
  for f in fmt.fields:
    if f.name == "GT": continue
    ctx.set_format_field(f, fmt, ints, floats)

  var alts = variant.format.genotypes(ints).alts
  for trio in ctx.trios:
      trio.fill("alts", alts, 1)

  for trio in ctx.trios:
    trio[0].duk.alias("kid")
    trio[1].duk.alias("dad")
    trio[2].duk.alias("mom")
    for i, dukex in ctx.expressions:
      if dukex.check:
        result.mgetOrPut(ctx.names[i], newSeq[string]()).add(samples[trio[0].i])


proc main*(dropfirst:bool=false) =
  let doc = """
variexpr -- variant expression for great good

Usage: variexpr [options] <expressions>...

Arguments:

    <expressions>...  as many name:expression pairs as required. an example would be:
    "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 \
                          && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10"
    this will be evaluated for every trio with kid, mom, dad set appropriately.
    other examples:

    "high_impact:/HIGH/.test(INFO.CSQ)"

    "rare_transmitted:(kid.alts > 0) && (dad.alts > 0 || mom.alts > 0) && kid.DP > 10 && mom.DP > 0 && INFO.AF < 0.01"

Options

  -v --vcf <path>       VCF/BCF
  -p --ped <path>       pedigree file with trio relations
  -o --out-vcf <path>   VCF/BCF
  --pass-only           only output variants that pass at least one of the filters [default: false]

  """

  var args: Table[string, docopt.Value]
  if dropfirst:
    var argv = commandLineParams()
    args = docopt(doc, argv=argv[1..argv.high])
  else:
    args = docopt(doc)

  if $args["--vcf"] == "nil":
    stderr.write_line "must specify the --vcf"
    quit doc
  if $args["--ped"] == "nil":
      stderr.write_line "must specify the --ped"
      quit doc
  if $args["--out-vcf"] == "nil":
    stderr.write_line "must specify the --out-vcf"
    quit doc
  var
    ivcf:VCF
    ovcf:VCF

  if not open(ivcf, $args["--vcf"], threads=3):
    quit "couldn't open:" & $args["--vcf"]

  var pass_only = bool(args["--pass-only"])

  var samples = parse_ped($args["--ped"])
  samples = samples.match(ivcf)
  var kids = newSeq[Sample]()

  for sample in samples:
    if sample.mom == nil or sample.dad == nil: continue
    kids.add(sample)

  stderr.write_line &"{kids.len} kids to be evaluated"

  if not open(ovcf, $args["--out-vcf"], mode="w"):
    quit "couldn't open:" & $args["--out-vcf"]

  ovcf.copy_header(ivcf.header)

  var tbl = initTable[string, string]()
  for e in @(args["<expressions>"]):
    var t = e.split(seps={':'}, maxsplit=1)
    if t.len != 2:
      quit "must specify name:expression pairs"
    tbl[t[0]] = t[1]
    if ovcf.header.add_info(t[0], ".", "String", &"added by variexpr with expression: '{t[1]}' from {$args[\"--vcf\"]}") != Status.OK:
      quit "error adding field to header"

  doAssert ovcf.write_header

  var ev = newEvaluator(kids, tbl)
  var t = cpuTime()
  var n = 10000
  var vcf_samples: seq[string] = ivcf.samples

  var i = 0
  for variant in ivcf:
    variant.vcf = ovcf
    i += 1
    if i mod n == 0:
      var secs = cpuTime() - t
      var persec = n.float64 / secs.float64
      stderr.write_line &"[variexpr] {i} {variant.CHROM}:{variant.start} evaluated {n} variants in {secs:.1f} seconds ({persec:.1f}/second)"
      t = cpuTime()
    var dns = ev.evaluate(variant, vcf_samples)
    if pass_only and dns.len == 0: continue
    for name, samples in dns:
      var ssamples = join(samples, ",")
      if variant.info.set(name, ssamples) != Status.OK:
        quit "error setting field:" & name

    doAssert ovcf.write_variant(variant)

  ovcf.close()
  ivcf.close()

when isMainModule:
  main()
