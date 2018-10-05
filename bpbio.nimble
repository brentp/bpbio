# Package

version       = "0.1.0"
author        = "Brent Pedersen"
description   = "collection of command-line tools for genomics"
license       = "MIT"
srcDir        = "src"
installExt    = @["nim"]
bin           = @["bpbio"]


# Dependencies

requires "nim >= 0.19.0", "plotly", "hts", "docopt"
