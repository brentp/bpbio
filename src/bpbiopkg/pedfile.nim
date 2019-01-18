import strutils
import math
import binaryheap
import sequtils
import strformat
import deques
import hashes
import sets
import tables
import hts/vcf

type
  Sample* = ref object
    family_id*: string
    id*: string
    paternal_id*: string
    maternal_id*: string
    mom*: Sample
    dad*: Sample
    sex*: int
    affected*: bool
    kids*: seq[Sample]
    i*:int

proc `$`*(s:Sample): string =
  return format(&"Sample(id:{s.id}, i:{s.i})")

iterator siblings*(s:Sample): Sample =
  ## report the samples with same mom or dad.
  ## TODO: look at both mom and dad.
  ## TODO: handle cases where mom and dad ids
  var kids = if s.mom != nil: s.mom.kids elif s.dad != nil: s.dad.kids else: newSeq[Sample]()
  for kid in kids:
    if kid != s: yield kid
  if s.mom == nil and s.dad == nil and s.maternal_id != "" and s.paternal_id != "":
    # TODO: send in seq[Samples] and find samples with same ids
    discard

proc spouse*(s:Sample): Sample {.inline.} =
  if s.kids.len == 0: return nil
  var k = s.kids[0]
  if k.dad == s: return k.mom
  return k.dad

proc addAncestors(q: var Deque[Sample], o:Sample) =
  if o.mom != nil:
    q.addLast(o.mom)
  if o.dad != nil:
    q.addLast(o.dad)
  if o.mom != nil:
    q.addAncestors(o.mom)
  if o.dad != nil:
    q.addAncestors(o.dad)

proc successors(o:Sample, samples: var seq[Sample]) =
  ## add kids and kids of kids
  for kid in o.kids:
    samples.add(kid)
  for kid in o.kids:
    successors(kid, samples)

proc successors(o:Sample): seq[Sample] =
  successors(o, result)

proc ancestors(o:Sample, samples:var seq[Sample]) =
  if o.mom != nil:
    samples.add(o.mom)
  if o.dad != nil:
    samples.add(o.dad)
  if o.mom != nil:
    ancestors(o.mom, samples)
  if o.dad != nil:
    ancestors(o.dad, samples)

proc ancestors(o:Sample): seq[Sample] =
  ancestors(o, result)

proc hash*(s:Sample): Hash =
  var h: Hash = 0

  h = h !& hash(s.family_id)
  h = h !& hash(s.id)
  result = !$h

proc lowest_common_ancestors(samples:seq[Sample], a:Sample, b:Sample): seq[Sample] =
  var all_ancestors = initSet[Sample](8)
  var common_ancestors = initSet[Sample](2)
  for i, target in @[a, b]:
    var parents = initSet[Sample](8)
    var q = initDeque[Sample]()

    q.addFirst(target)
    while q.len > 0:
      var n = q.popFirst
      if parents.contains(n): continue
      # A predecessor of n is a node m such that there exists a directed edge from m to n.
      # A successor of n is a node m such that there exists a directed edge from n to m.
      #
      # TODO: might need to flip addKids and ancestors.
      q.addAncestors(n)
      parents.incl(n)

    if i == 0:
      common_ancestors.incl(parents)
    else:
      common_ancestors = common_ancestors.intersection(parents)
    all_ancestors.incl(parents)


  for p in common_ancestors:
    for child in p.successors:
      if not common_ancestors.contains(child) and child in all_ancestors:
        result.add(p)
        break

type dsample = object
  s: Sample
  dist: int

proc dijkstra(samples:seq[Sample], a:Sample, b:Sample): int =
  ## return the shortest path from a to b. in this case a must be
  ## an predecessor (older than) b.
  ## see: https://github.com/mburst/dijkstras-algorithm/blob/master/dijkstras.py
  var dsamples = newSeqOfCap[dsample](4)
  var previous = initTable[Sample, Sample](4)
  if a.family_id != b.family_id: return -1

  var h = newHeap[dsample]() do (a, b: dsample) -> int:
    return a.dist - b.dist

  for v in samples:
    if v.family_id != a.family_id: continue
    if v == a: dsamples.add(dsample(s:v, dist:0))
    else: dsamples.add(dsample(s:v, dist:samples.len + 2))
    h.push(dsamples[dsamples.high])

  while h.size > 0:
    var dsmallest = h.pop()
    var smallest = dsmallest.s
    if smallest == b:
      if not previous.contains(smallest): return -1
      while previous.contains(smallest):
        result += 1
        smallest = previous[smallest]
      return

    if dsmallest.dist == samples.len + 2: break

    for okid in dsmallest.s.kids:
      # TODO: make this more efficient.
      var kid:dsample
      for s in dsamples:
        if s.s == okid:
          kid = s
          break
      var alt = dsmallest.dist + 1
      if alt < kid.dist:
        # TODO: make sure kid with updated dist gets into the heap
        kid.dist = alt
        previous[okid] = dsmallest.s

        var heapSamples = toSeq(h.items)

        h = newHeap[dsample](proc(a, b: dsample):  int =
          return a.dist - b.dist
        )
        for s in heapSamples:
          if s.s.id == okid.id:
            h.push(kid)
          else:
            h.push(s)
  return -1


proc relatedness*(a:Sample, b:Sample, samples: seq[Sample]): float64 =
  ## coefficient of relatedness of 2 samples given their family.
  ## to differentiates siblings from parent-child relationships,
  ## siblings are return with a value of 0.49
  ## samples from different families are given a value of -1.
  if a.family_id != b.family_id: return -1
  if a.dad == b or a.mom == b: return 0.5
  if b.dad == a or b.mom == a: return 0.5

  if a notin samples: return -1'f64
  if b notin samples: return -1'f64

  var lca = lowest_common_ancestors(samples, b, a)
  if lca.len == 0: return 0
  var
    amax: int = 1000 # hacky use of 1000 and empty value.
    bmax: int = 1000
    n = 0

  for anc in lca:
    if anc != a:
      var d = dijkstra(samples, anc, a)
      if d > 0:
        amax = min(amax, d)
    if anc != b:
      var d = dijkstra(samples, anc, b)
      if d > 0:
        bmax = min(bmax, dijkstra(samples, anc, b))

  # if one of the samples is in the set, then their distance is 1.
  if a in lca:
    amax = 1
  if b in lca:
    bmax = 1

  if amax != 1000: n += 1 else: amax = 0
  if bmax != 1000: n += 1 else: bmax = 0

  # subtract 2 because the path includes the end sample which we don't need.
  return n.float64 * pow(2'f64, -float64(amax + bmax))

proc parse_ped*(path: string, verbose:bool=true): seq[Sample] =
  result = new_seq_of_cap[Sample](10)

  var look = newTable[string,Sample]()

  for line in lines(path):
    if line.len > 0 and line[0] == '#': continue
    if line.strip().len == 0: continue
    var toks = line.strip().split('\t')
    if toks.len < 6:
      stderr.write_line "[pedfile] error: expected at least 5 tab-delimited columns in ped file: " & path
      stderr.write_line "[pedfile] error: line was:" & $toks


    var s = Sample(family_id: toks[0], id: toks[1], kids:new_seq[Sample](), paternal_id: toks[2], maternal_id:toks[3], i: -1)
    s.affected = toks[5] == "2"
    s.sex = if toks[4] == ".": -9 else: parseInt(toks[4])
    result.add(s)
    look[s.id] = s

  for s in result:
    if s.paternal_id in look:
      s.dad = look[s.paternal_id]
      s.dad.kids.add(s)
    elif verbose and not (s.paternal_id in @[".", "-9", "", "0"]):
      stderr.write_line &"[pedfile] paternal_id: \"{s.paternal_id}\" referenced for sample {s.id} not found"
    if s.maternal_id in look:
      s.mom = look[s.maternal_id]
      s.mom.kids.add(s)
    elif verbose and not (s.maternal_id in @[".", "-9", "", "0"]):
      stderr.write_line &"[pedfile] maternal_id: \"{s.maternal_id}\" referenced for sample {s.id} not found"

proc match*(samples: seq[Sample], vcf:var VCF, verbose:bool=true): seq[Sample] =
  ## adjust the VCF samples and the samples to match
  var si = newTable[string,int]()
  for i, s in vcf.samples:
    si[s] = i
  # get samples in same order as VCF
  result = new_seqOfCap[Sample](len(si))
  var samplesById = newTable[string,Sample]()
  for sample in samples:
    sample.i = -1
    if not (sample.id in si):
      if verbose:
        stderr.write_line(&"[pedfile] {sample.id} not found in VCF samples")
      continue

    sample.i = si[sample.id]
    samplesById[sample.id] = sample

  var sampleList = newSeq[string]()
  for s in vcf.samples:
    if not (s in samplesById):
      if verbose:
        stderr.write_line(&"[pedfile] {s} from VCF not found in samples")
      continue
    sampleList.add(s)

  vcf.set_samples(sampleList)
  for i, s in vcf.samples:
    var sample = samplesById[s]
    sample.i = i
    result.add(sample)

when isMainModule:
  import algorithm
  import unittest
  import times

  var samples = parse_ped("tests/testa.ped")
  var ovcf:VCF
  if not open(ovcf, "tests/test.vcf"):
    quit "bad"

  suite "pedfile test suite":

    test "that samples match those in vcf":
      var osamples = samples.match(ovcf)
      for i, s in osamples:
        check s.i == i
        check s.id == ovcf.samples[i]

    test "that reversed samples match those in vcf":
      var osamples = samples.reversed
      osamples = osamples.match(ovcf)
      for i, s in osamples:
        check s.i == i
        check s.id == ovcf.samples[i]

    test "siblings":
      check samples[0].id == "101976"
      for sib in samples[0].siblings:
        check sib.dad == samples[0].dad
        check sib.mom == samples[0].mom


    test "lowest common ancestor":
      var k1 = Sample(family_id:"1", id:"kid1")
      var k2 = Sample(family_id:"1", id:"kid2")
      var mom = Sample(family_id:"1", id:"mom")
      var dad = Sample(family_id:"1", id:"dad")
      k1.mom = mom
      k2.mom = mom
      k1.dad = dad
      k2.dad = dad
      mom.kids.add(@[k1, k2])
      dad.kids.add(@[k1, k2])

      echo lowest_common_ancestors(@[k1, k2, mom, dad], k1, k2)

    test "dijkstra and relatedness":
      var k1 = Sample(family_id:"1", id:"kid1")
      var k2 = Sample(family_id:"1", id:"kid2")
      var mom = Sample(family_id:"1", id:"mom")
      var dad = Sample(family_id:"1", id:"dad")
      var uncle = Sample(family_id: "1", id:"uncle")
      var cousin = Sample(family_id: "1", id:"cousin")
      var gma = Sample(family_id:"1", id:"gma")
      var ggma = Sample(family_id:"1", id:"ggma")
      var unrel = Sample(family_id:"1", id:"un")
      var extern = Sample(family_id:"xxx", id:"extern")
      k1.mom = mom
      k2.mom = mom
      k1.dad = dad
      k2.dad = dad
      dad.mom = gma
      uncle.mom = gma
      uncle.kids.add(cousin)
      cousin.dad = uncle
      gma.kids.add(@[dad, uncle])
      mom.kids.add(@[k1, k2])
      dad.kids.add(@[k1, k2])
      ggma.kids.add(gma)
      gma.mom = ggma
      var fam = @[k1, k2, mom, dad, gma, ggma, cousin, uncle, unrel]
      check 2 == dijkstra(fam, gma, k1)
      check 1 == dijkstra(fam, mom, k1)
      check 3 == dijkstra(fam, ggma, k1)
      check 1 == dijkstra(fam, ggma, gma)
      check 2 == dijkstra(fam, ggma, dad)
      check -1 == dijkstra(fam, ggma, mom)

      check relatedness(uncle, k1, fam) == 0.25
      check relatedness(dad, k1, fam) == 0.5
      check relatedness(k1, k2, fam) == 0.5
      check relatedness(k2, mom, fam) == 0.5
      check relatedness(k2, gma, fam) == 0.25
      check relatedness(k2, ggma, fam) == 0.125
      check relatedness(extern, gma, fam) == -1'f64
      check relatedness(unrel, gma, fam) == 0.0'f64
      check relatedness(k1, cousin, fam) == 0.125

    test "relatedness":
      var k1 = Sample(family_id:"1", id:"kid1")
      var dad = Sample(family_id:"1", id:"kid1")
      k1.dad = dad
      dad.kids.add(k1)
      var fam = @[k1, dad]

      check 0.5 == relatedness(k1, dad, fam)

    test "relatedness with no relations":
      var a = Sample(family_id:"1", id:"a")
      var b = Sample(family_id:"2", id:"b")
      check relatedness(a, b, @[a, b]) == -1'f64
      b.family_id = "1"
      check relatedness(a, b, @[a, b]) == -0'f64


    test "relatedness with 603 ceph samples":
      var samples = parse_ped("tests/ceph.ped")
      var t = cpuTime()
      var n = 0

      for i, sampleA in samples[0..<samples.high]:
        for j, sampleB in samples[i + 1..samples.high]:
          var rel = relatedness(sampleA, sampleB, samples)
          doAssert -1'f64 <= rel
          if rel > 0.5:
            echo &"{sampleA} {sampleB} {rel}"
          n += 1
          if n mod 50000 == 0: echo "tested:", n

      echo &"time for {n} calculations: {cpuTime() - t:.1f}"
