import strutils
import strformat
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

proc is_grandkid*(s:Sample): bool {.inline.} =
  return s.mom != nil and s.mom.dad != nil and s.dad != nil and s.mom.mom != nil and s.dad.dad != nil and s.dad.mom != nil

proc is_f1*(s:Sample): bool {.inline.} =
  return s.mom != nil and s.dad != nil and s.kids.len > 0

proc spouse*(s:Sample): Sample {.inline.} =
  if s.kids.len == 0: return nil
  var k = s.kids[0]
  if k.dad == s: return k.mom
  return k.dad

proc parse_ped*(path: string, verbose:bool=true): seq[Sample] =
  result = new_seq_of_cap[Sample](10)

  var look = newTable[string,Sample]()

  for line in lines(path):
    if line[0] == '#': continue
    var toks = line.strip().split('\t')

    var s = Sample(family_id: toks[0], id: toks[1], kids:new_seq[Sample](), paternal_id: toks[2], maternal_id:toks[3], i: -1)
    s.affected = toks[5] == "2"
    s.sex = if toks[4] == ".": -9 else: parseInt(toks[4])
    result.add(s)
    look[s.id] = s

  for s in result:
    if s.paternal_id in look:
      s.dad = look[s.paternal_id]
      s.dad.kids.add(s)
    elif verbose and not (s.paternal_id in @[".", "-9", ""]):
      stderr.write_line &"[pedfile] paternal_id: \"{s.paternal_id}\" referenced for sample {s.id} not found"
    if s.maternal_id in look:
      s.mom = look[s.maternal_id]
      s.mom.kids.add(s)
    elif verbose and not (s.maternal_id in @[".", "-9", ""]):
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
