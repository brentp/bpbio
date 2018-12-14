import ./bpbiopkg/variexpr
import os
import unittest

suite "hwe scoring":
  test "hwe is ok":
    var counts = [30, 55, 15, 0]
    check counts.aaf == 0.425
    check abs(counts.hwe_score(counts.aaf) - 1.57) < 0.01

