# This is just an example to get you started. A typical hybrid package
# uses this file as the main entry point of the application.

import bpbiopkg/plotsvvcf
from bpbiopkg/info import version
import strformat
import tables
import os

type pair = object
    f: proc()
    description: string

var dispatcher = {
  "plot-sv-vcf": pair(f:plotsvvcf.main, description:"make a plot of SV types across samples for a multi-sample SV VCF")
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

  dispatcher[args[0]].f()
