import plotly
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

proc main*() =
  let doc = format("""

    Usage: bpbio plot-sv-vcf <vcf>

Arguments:

    <vcf>    a structural variant VCF

    """)
  let args = docopt(doc)
  var vcf:VCF

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

  for v in smalls.mvalues:
    v.text = toText(vcf.samples)
  for v in larges.mvalues:
    v.text = toText(vcf.samples)


  var svtype: string
  var gts = newSeq[int32](n)
  var tr: Trace[int]
  var k = 0
  for v in vcf:
    k += 1
    if k mod 10000 == 0:
      echo $k & " > " & $v.CHROM & ":" & $v.start
      #echo larges["DUP"].ys
    if v.info.get("SVTYPE", svtype) != Status.OK:
       quit "no svtype for:" & $v
    var svlen = v.stop - v.start
    var nalts:seq[int8] = v.format.genotypes(gts).alts

    if svlen > 1000:
      tr = larges[svtype]
    else:
      tr = smalls[svtype]
    for i, na in nalts:
      if na == 0 or na == -1:
        continue
      tr.ys[i] += 1

  var tmpl_hdr = """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>$title</title>
   <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    """
  var tmpl_body = """
      <div id="plot$i"></div>
      <script>
      Plotly.newPlot('plot$i', $data, $layout)
      </script>
    """

  var html = tmpl_hdr % ["title", $args["<vcf>"]]

  for i, pt in @["DEL", "DUP", "INV", "BND"]:
    var layout = Layout(width: 1200, height: 400,
                    yaxis: Axis(title:pt),
                    xaxis: Axis(title: "sample", hideticklabels:true),
                    barmode: BarMode.Stack,
                    autosize: false)
    var jsons = mapIt(@[smalls[pt], larges[pt]], it.json(as_pretty = false))
    var j = "[" & join(jsons, ",") & "]"

    var h = tmpl_body % ["i", $i, "data", j, "layout", $(%layout)]
    html &= h

    #$Plot[int](layout: layout, traces: @[smalls[pt], larges[pt]]).show()
  html &= "</body></html>"
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
