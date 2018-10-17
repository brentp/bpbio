import unittest

import "bpbiopkg/homsv"

suite "flip suite":
    test "that flip works":
        var a = region_tid(achrom: 1, astart:22, astop:33, bchrom: 2, bstart:222, bstop:555)

        var b = flip(a)
        check b.achrom == 2
        check b.astart == 222
        check b.astop == 555
        check b.bchrom == 1
        check b.bstart == 22
        check b.bstop == 33
