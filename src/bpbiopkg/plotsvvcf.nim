import plotly
import chroma
import browsers
import os
import hts
import docopt
import json
import strutils
import sequtils

proc toText(samples:seq[string]): seq[string] =
  result = newSeq[string](samples.len)
  for i, s in samples:
    result[i] = "sample:" & $s

var afsColor = Color(r: 0.5, g: 0.5, b: 0.5, a:1.0)

proc get_bnd_mate_pos*(a:string, vchrom:string): int {.inline.} =
    if not (':' in a): return -1
    var tmp = a.split(':')
    var left = tmp[0]
    var i = 0
    while i < 3:
      if left[i] in {'[', ']'}:
        break
      i += 1

    var chrom = left[i+1..left.high]
    if chrom != vchrom: return -1

    var right = tmp[1]
    i = 0
    while right[i].isdigit:
        i += 1
    result = parseInt(right[0..<i])

proc get_bnd_mate_pos(variant:Variant): int {.inline.} =
    return get_bnd_mate_pos(variant.ALT[0], $variant.CHROM)

proc main*() =
  let doc = format("""

    Usage: bpbio plot-sv-vcf [options] <vcf>

Arguments:

   <vcf>    a structural variant VCF

Options:

  -t --threads <int>       number of BAM decompression threads [default: 2]
  -s --size-cutoff <int>   size below which events are considered small [default: 300]


    """)
  let args = docopt(doc)
  var vcf:VCF
  let threads = parseInt($args["--threads"])
  let size = parseInt($args["--size-cutoff"])

  if not open(vcf, $args["<vcf>"], threads=2):
     quit "couldn't open vcf:" & $args["<vcf>"]

  var n = vcf.n_samples

  var
    smalls = {
       "DEL": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"small deletions", ys: newSeq[int](n)),
       "DUP": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"small duplications", ys: newSeq[int](n)),
       "INV": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"small inversion", ys: newSeq[int](n)),
       "BND": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"small BNDs", ys: newSeq[int](n)),
    }.toTable
    larges = {
       "DEL": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"large deletions", ys: newSeq[int](n)),
       "DUP": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"large duplications", ys: newSeq[int](n)),
       "INV": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"large inversion", ys: newSeq[int](n)),
       "BND": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"large BNDs", ys: newSeq[int](n)),
    }.toTable
    inters = {
       "BND": Trace[int](`type`: PlotType.Bar, opacity: 0.8, name:"interchromosomal BNDs", ys: newSeq[int](n)),
    }.toTable
    type_counts = {
       "DEL": 0,
       "DUP": 0,
       "INV": 0,
       "BND": 0,
    }.toTable
    AFS_counts = {
      "DEL": Trace[int](`type`: PlotType.Histogram, xs: newSeqOfCap[int](4096), marker:Marker[int](color:(@[afsColor]))),
      "DUP": Trace[int](`type`: PlotType.Histogram, xs: newSeqOfCap[int](4096), marker:Marker[int](color:(@[afsColor]))),
      "INV": Trace[int](`type`: PlotType.Histogram, xs: newSeqOfCap[int](4096), marker:Marker[int](color:(@[afsColor]))),
      "BND": Trace[int](`type`: PlotType.Histogram, xs: newSeqOfCap[int](4096), marker:Marker[int](color:(@[afsColor]))),
    }.toTable

  for v in smalls.mvalues:
    v.text = toText(vcf.samples)
  for v in larges.mvalues:
    v.text = toText(vcf.samples)
  for v in inters.mvalues:
    v.text = toText(vcf.samples)


  var svtype: string
  var gts = newSeq[int32](n)
  var tr: Trace[int]
  var k = 0
  for v in vcf:
    k += 1
    if k mod 50_000 == 0:
      echo $k & " > " & $v.CHROM & ":" & $v.start
      #echo larges["DUP"].ys
    if v.info.get("SVTYPE", svtype) != Status.OK:
       quit "no svtype for:" & $v
    var svlen = v.stop - v.start
    var nalts:seq[int8] = v.format.genotypes(gts).alts
    type_counts[svtype].inc
    if svtype == "BND":
      var p = get_bnd_mate_pos(v)
      if p == -1:
        svlen = -1
      else:
        svlen = (v.start - p).abs
    if svlen == -1:
      tr = inters[svtype]
    elif svlen > size:
      tr = larges[svtype]
    else:
      tr = smalls[svtype]
    var n_samples_with_event = 0
    for i, na in nalts:
      if na == 0 or na == -1:
        continue
      n_samples_with_event += 1
      tr.ys[i] += 1
    AFS_counts[svtype].xs.add(n_samples_with_event)

  var tmpl_hdr = """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>$title</title>
   <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
<table>
    """
  var tmpl_body = """
<td>
      <div id="plot$i"></div>
      <script>
      var pdata_$i = $data;
      var playout_$i = $layout;
      playout_$i.legend = {orientation: 'h', x:0, y:0};
      if(Number($i) >= 10) {
        playout_$i.yaxis.log = true;
        playout_$i.yaxis.type = "log";
      }
      var plot$i = Plotly.newPlot('plot$i', pdata_$i, playout_$i)
      </script>
</td>
    """

  var html = tmpl_hdr % ["title", $args["<vcf>"]]

  for i, pt in @["DEL", "DUP", "INV", "BND"]:
    var layout = Layout(width: 950, height: 350,
                    yaxis: Axis(title:pt & " (" & $type_counts[pt] & ")"),
                    xaxis: Axis(title: "sample", hideticklabels:true),
                    barmode: BarMode.Stack,
                    autosize: false)
    var jsons:seq[string]
    if pt in inters:
      jsons = mapIt(@[smalls[pt], larges[pt], inters[pt]], it.json(as_pretty = false))
    else:
      jsons = mapIt(@[smalls[pt], larges[pt]], it.json(as_pretty = false))
    var j = "[" & join(jsons, ",") & "]"

    var h = tmpl_body % ["i", $i, "data", j, "layout", $(%layout)]
    html &= "<tr>" & h

    # AFS plot
    layout = Layout(width: 950, height: 350,
                    xaxis: Axis(title: "# of samples with " & pt & "s"),
                    yaxis: Axis(title: "count of variants"),
                    autosize: false)
    jsons = @[AFS_counts[pt].json(as_pretty=false)]
    j = "[" & join(jsons, ",") & "]"

    h = tmpl_body % ["i", $(i + 10), "data", j, "layout", $(%layout)]
    html &= h & "</tr>"


  html &= "</table></body></html>"
  var
     f: File
  if not open(f, "svvcf.html", fmWrite):
    quit "could not open file for json"
  f.write(html)
  f.close()
  sleep(100)
  openDefaultBrowser("svvcf.html")


when isMainModule:
  main()
