import docopt
import hts
import genoiser
import strutils


proc ifun(aln:Record, posns: var seq[mrange]) =
  if aln.mapping_quality == 0:
    return
  var f = aln.flag
  if f.unmapped or f.secondary or f.qcfail or f.dup: return
  posns.add((aln.start, aln.stop, 1))


proc find(targets: seq[Target], chrom: string): Target =
    for t in targets:
        if t.name == chrom: return t
    return nil

type cov[T] = tuple[sum:seq[T], variance:seq[float64]]

proc fastcov*[T](paths: seq[string], fasta:string, threads:int, chrom: string): cov[T] =
  ## fast coverage (doesn't do cigar operations or correct for pair overlap) for a single
  ## chromosome. can be used for multiple BAM/CRAM files to get the sum of coverage.

  var depths = Fun[T](f:ifun)
  var single = Fun[T](f:ifun)
  var S: seq[float64]
  var M: seq[float64]
  for i, path in paths:
    var bam:Bam
    open(bam, path, threads=threads, fai=fasta, index=true)
    if bam == nil:
      quit "couldn't open bam:" & path
    discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 511)
    var target = bam.hdr.targets.find(chrom)
    if target == nil:
        quit "couldn't find " & chrom & "in " & path
    if i == 0:
        depths.values = newSeq[T](target.length.int + 1)
        S = newSeq[float64](target.length.int + 1)
        M = newSeq[float64](target.length.int + 1)
    single.values = newSeq[T](target.length.int + 1)
    discard genoiser[T](bam, @[single], target.name, 0, target.length.int, cumsum=false)
    for k, v in single.values:
        if depths.values[k].int64 + v.int64 < T.high:
          depths.values[k] += v
        else:
          depths.values[k] = T.high
        if i == 0:
          M[k] = v.float64
          S[k] = 0'f64
        else:
          var oldm = M[k]
          var olds = S[k]
          M[k] = oldm + (v.float64 - oldm) / (i + 1).float64
          S[k] = olds + (v.float64 - oldm)*(v.float64 - M[k])

    bam.close()

  depths.values.cumulative_sum()
  return (depths.values, S)


proc main*() =
  let doc = format("""
Very fast coverage calculation.

    Usage: bpbio fastcov [options] <fasta> <BAM>...

Arguments:

   <BAM>    indexed BAM or CRAM files

Options:

  -t --threads <int>       number of BAM decompression threads [default: 2]
  -s --size-cutoff <int>   cigar events less than this size are ignored [default: 20]
  -c --chrom <string>      chromosome to limit depth [default: 1]
    """)
  let args = docopt(doc)
  var bam:BAM
  let threads = parseInt($args["--threads"])
  let size = parseInt($args["--size-cutoff"])
  let fasta = $args["<fasta>"]

  var dp = fastcov[int32](@(args["<BAM>"]), fasta, threads, $args["--chrom"])

when isMainModule:
  main()
