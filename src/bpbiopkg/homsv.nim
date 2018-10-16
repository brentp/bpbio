import docopt
import algorithm
import strformat
import hts
import genoiser
import strutils
import tables
import duphold

type Strand* {.pure, size:1.} = enum
  Minus = 0
  Plus = 1
  PlusMinus = 2
  PlusPlus = 3
  MinusPlus = 4
  MinusMinus = 5

type SignalType* {.pure, size:1.} = enum
  Splitter = 0
  Pair = 1
  CigarEvent = 2

# 16 bytes
type Signal = object
  left_tid: uint16
  right_tid: uint16
  left: uint32
  right: uint32
  typ: SignalType
  strands: Strand
  value: uint16 # can store, e.g. mismatches and quality


type region_t* = ref object
    achrom*: string
    astart*: uint32
    astop*: uint32
    bchrom*: string
    bstart*:uint32
    bstop*: uint32
    percent_id*: float32

type region_tid* = ref object
    achrom*: uint16
    bchrom*: uint16
    astart*: uint32
    astop*: uint32
    bstart*:uint32
    bstop*: uint32
    percent_id*: float32

proc `$`(r:region_tid): string =
    return &"{r.achrom}:{r.astart}-{r.astop} ... {r.bchrom}:{r.bstart}-{r.bstop}"

proc chain_line_to_region(line: string, stripchr:bool): region_t =
  # 585     2597464 chr1    249250621       10000   39170   chr12   133851895       -       133756271       133787376       4031    89.7
  var
    cse = line.strip().split('\t', 15)
  result = region_t(achrom: cse[2], bchrom: cse[6])
  result.astart = parse_int(cse[4]).uint32
  result.astop = parse_int(cse[5]).uint32
  result.bstart = parse_int(cse[9]).uint32
  result.bstop = parse_int(cse[10]).uint32
  result.percent_id = parse_float(cse[12]).float32
  if stripchr:
    if result.achrom.startswith("chr"):
      result.achrom = result.achrom[3..result.achrom.high]
    if result.bchrom.startswith("chr"):
      result.bchrom = result.bchrom[3..result.bchrom.high]
  if result.achrom != result.bchrom: return
  if result.astart > result.bstart:
    var t = result.astart
    result.astart = result.bstart
    result.bstart = t
    t = result.astop
    result.astop = result.bstop
    result.bstop = t
  #if result.astop > result.bstart:
  #    quit "bad region:" & $result[]

proc flip*(a: region_tid): region_tid =
    result = region_tid(achrom: a.bchrom, astart:a.bstart, astop:a.bstop,
                      bchrom: a.achrom, bstart:a.astart, bstop:a.astop, percent_id:a.percent_id)

proc toTid*(r: region_t, tids: seq[string]): region_tid =
    var atid = tids.find(r.achrom)
    if atid == -1: return nil
    var btid = tids.find(r.bchrom)
    if btid == -1: return nil
    result = region_tid(achrom: atid.uint16, astart:r.astart, astop:r.astop,
                        bchrom: btid.uint16, bstart:r.bstart, bstop:r.bstop,
                        percent_id:r.percent_id)


proc bed_to_table(bed: string, stripchr:bool, tids:seq[string]): TableRef[uint16, seq[region_tid]] =
  result = newTable[uint16, seq[region_tid]]()
  var kstr = kstring_t(l:0, m: 0, s: nil)
  var hf = hts_open(cstring(bed), "r")
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if kstr.s[0] == 't' and ($kstr.s).startswith("track "):
      continue
    if kstr.s[0] == '#':
      continue
    var v = chain_line_to_region($kstr.s, stripchr).toTid(tids)
    if v == nil: continue

    discard result.hasKeyOrPut(v.achrom, new_seq[region_tid]())
    result[v.achrom].add(v)
    if v.achrom != v.bchrom:
        var o = v.flip
        discard result.hasKeyOrPut(o.achrom, new_seq[region_tid]())
        result[o.achrom].add(o)

  # since it is read into mem, can also well sort.
  for chrom, ivs in result.mpairs:
      sort(ivs, proc (a, b: region_tid): int = int(a.astart) - int(b.astart))

  hts.free(kstr.s)

proc median_depth(depth:var seq[int16], start:uint32, stop:uint32, m:var MedianStats) {.inline.} =
    m.clear
    for i in start..stop:
        m.add(depth[i.int].int)

proc main*() =
  let doc = format("""
find SVs between self-similar regions.

    Usage: bpbio homsv [options] <fasta> <BAM>

Arguments:

   <fasta>  path to fasta file
   <BAM>    path to BAM or CRAM file

Options:

  -t --threads <int>       number of BAM decompression threads [default: 2]
  -c --self-chain <path>   path to selfchain.bed.gz
  -s --strip-chr           strip "chr" prefix from the self-chain file.
    """)
  let args = docopt(doc)
  var bam:BAM
  let
    threads = parseInt($args["--threads"])
    fasta = $args["<fasta>"]
    strip = args["--strip-chr"] == true
  if $args["--self-chain"] == "nil":
      quit "--self-chain argument required"

  open(bam, $args["<BAM>"], threads=threads, fai=fasta, index=true)
  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 511)

  var isize = bam.get_isize_distribution()
  var i98 = isize.percentile(98)

  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 511 + SAM_AUX.int)

  var tids = newSeq[string]()
  for t in bam.hdr.targets:
      tids.add(t.name)

  var regions = bed_to_table($args["--self-chain"], strip, tids)

  for t in bam.hdr.targets:
    if not (t.tid.uint16 in regions):
      if t.length.int > 100000:
        stderr.write_line &"{t.name} not found in regions"
      continue

    var signals = newSeqOfCap[Signal](32768)

    proc ifun(aln:Record, posns: var seq[mrange]) =
      if aln.mapping_quality == 0:
        return
      var f = aln.flag
      if f.unmapped or f.secondary or f.qcfail or f.dup or f.supplementary: return
      posns.add((aln.start, aln.stop, 1))
      if aln.mate_pos < aln.start: return
      # TODO: handle interchromosomals.
      if aln.tid != aln.mate_tid: return
      #if aln.isize < 0: return
      #if not f.read1: return
      if aln.isize.abs < i98 and f.reverse != f.mate_reverse:
          return

      var s = Signal(left_tid: aln.tid.uint16, right_tid: aln.tid.uint16)
      s.left = min(aln.stop, aln.mate_pos).uint32
      s.right = max(aln.stop, aln.mate_pos).uint32
      s.typ = SignalType.Pair

      if aln.isize >= i98:
          return
      if f.reverse and not f.mate_reverse:
        s.strands = Strand.MinusPlus
      elif f.reverse and f.mate_reverse:
        s.strands = Strand.MinusMinus
      elif not (f.reverse or f.mate_reverse):
        s.strands = Strand.PlusPlus
      else:
          quit "unhandled alingment:" & aln.tostring
      #echo f.int
      signals.add(s)

    var depths = Fun[int16](f:ifun, values: newSeq[int16](t.length.int + 1))

    if not genoiser[int16](bam, @[depths], t.name, 0, t.length.int):
        continue
    stderr.write_line &"{t.name} found in regions"

    var d = depths.values
    var m = MedianStats()
    d.median_depth(0.uint32, d.high.uint32, m)
    var chrom_med = m.median(skip_zeros=false)

    for r in regions[t.tid.uint16]:
      if r.achrom != r.bchrom: continue
      d.median_depth(r.astop, r.bstart, m)
      var lm = m.median(skip_zeros=false)
      echo &"local median: {lm} global: {chrom_med}, relative depth: { lm.float64 / chrom_med.float64} region: {r}"


    stderr.write_line signals.len
    return

when isMainModule:
  main()
