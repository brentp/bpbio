# This is just an example to get you started. A typical hybrid package
# uses this file as the main entry point of the application.

import bpbiopkg/plotsvvcf
from bpbiopkg/info import version
#import bpbiopkg/fastcov
#import bpbiopkg/homsv
#import bpbiopkg/variexpr
#export fastcov
import strformat
import tables
import os

type pair = object
    f: proc(dropfirst:bool)
    description: string

var dispatcher = {
  "plot-sv-vcf": pair(f:plotsvvcf.main, description:"make a plot of SV types across samples for a multi-sample SV VCF"),
  #"variexpr": pair(f:variexpr.main, description:"simple expression to filter or annotate VCFs"),
  #"homsv": pair(f:homsv.main, description:"look for depth changes in self-chains or homologous regions"),
  #"homsv-merge": pair(f:homsv.merge, description:"merge output from homsv"),
  }.toTable

when isMainModule:
  stderr.write_line "bpbio version:" & version
  var args = commandLineParams()

  if len(args) == 0 or not (args[0] in dispatcher):
    for k, v in dispatcher:
      echo &"{k}:   {v.description}"
    if len(args) > 0 and not (args[0] in dispatcher):
        echo &"unknown program '{args[0]}'"
    quit 1

  dispatcher[args[0]].f(true)
