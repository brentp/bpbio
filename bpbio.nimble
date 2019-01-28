# Package
version       = "0.1.1"
author        = "Brent Pedersen"
description   = "collection of command-line tools and small utility functions for genomics"
license       = "MIT"
srcDir        = "src"
installExt    = @["nim"]
bin           = @["bpbio"]

# Dependencies

requires "nim >= 0.19.0", "plotly", "hts", "docopt", "genoiser >= 0.2.3", "binaryheap"
