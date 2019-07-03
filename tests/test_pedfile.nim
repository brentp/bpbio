import bpbiopkg/pedfile
import unittest


var samples = parse_ped("tests/distant.ped")

proc get_by_id(samples:seq[Sample], sid:string): Sample =
  for s in samples:
    if s.id == sid: return s
  raise newException(KeyError, sid)

#8125 8130
#
suite "extended pedigree":
    test "double first cousins":
      echo "8141 .. 8130"
      var a = samples.get_by_id("8141")
      var b = samples.get_by_id("8130")
      #var a = samples.get_by_id("8109")
      #var b = samples.get_by_id("8110")

      check 0.25 == a.relatedness(b)

    test "uncle":
      var a = samples.get_by_id("8135")
      var b = samples.get_by_id("8141")

      check 0.25 == a.relatedness(b)

    test "sibs":

      var kida = Sample(family_id:"1", id:"kida")
      var kidb = Sample(family_id:"1", id:"kidb")
      var dad = Sample(family_id:"1", id:"dad")
      var mom = Sample(family_id:"1", id:"mom")
      var uncle = Sample(family_id:"1", id:"uncle")

      var gma = Sample(family_id:"1", id:"gma")
      var gpa = Sample(family_id:"1", id:"gpa")

      var ggma = Sample(family_id:"1", id:"ggma")
      var ggpa = Sample(family_id:"1", id:"ggpa")

      dad.kids = @[kida, kidb]
      mom.kids = @[kida, kidb]
      kida.dad = dad
      kidb.dad = dad
      kida.mom = mom
      kidb.mom = mom

      dad.mom = gma
      dad.dad = gpa
      uncle.mom = gma
      uncle.dad = gpa

      gma.mom = ggma
      gma.dad = ggpa
      ####

      gma.kids = @[dad, uncle]
      gpa.kids = @[dad, uncle]

      var samples = @[kida, kidb, dad, mom, gma, gpa]

      check 0.49 == kida.relatedness(kidb)
      check 0.25 == kida.relatedness(uncle)

      # grandparents are full sibs
      gpa.mom = ggma
      gpa.dad = ggpa

      check 0.5 < kida.relatedness(kidb)

    test "ext pedigree":
      var samples = parse_ped("tests/test-ext.ped")
      var a = samples.get_by_id("8467")
      var b = samples.get_by_id("8472")
      #echo samples
      check 0.5 == a.relatedness(b)

    test "sibs":
      var k1 = Sample(family_id:"1", id:"kid1")
      var k2 = Sample(family_id:"1", id:"kid2")
      var mom = Sample(family_id:"1", id:"mom")
      var dad = Sample(family_id:"1", id:"dad")
      var uncle = Sample(family_id: "1", id:"uncle")
      var gma = Sample(family_id: "1", id:"gma")
      var gpa = Sample(family_id: "1", id:"gpa")

      k1.mom = mom
      k2.mom = mom
      k1.dad = dad
      k2.dad = dad

      mom.mom = gma
      mom.dad = gpa
      uncle.mom = gma
      uncle.dad = gpa
      check k1.relatedness(k2) == 0.49
      check k1.relatedness(uncle) == 0.25
      check uncle.relatedness(k1) == 0.25
      check dad.relatedness(k1) == 0.5

    test "relatedness":
      var k1 = Sample(family_id:"1", id:"kid1")
      var k2 = Sample(family_id:"1", id:"kid2")
      var mom = Sample(family_id:"1", id:"mom")
      var dad = Sample(family_id:"1", id:"dad")
      var uncle = Sample(family_id: "1", id:"uncle")
      var cousin = Sample(family_id: "1", id:"cousin")
      var gma = Sample(family_id:"1", id:"gma")
      var gpa = Sample(family_id:"1", id:"gma")
      var ggma = Sample(family_id:"1", id:"ggma")
      var ggpa = Sample(family_id:"1", id:"ggpa")
      var unrel = Sample(family_id:"1", id:"un")
      var extern = Sample(family_id:"xxx", id:"extern")
      k1.mom = mom
      k2.mom = mom
      k1.dad = dad
      k2.dad = dad
      dad.mom = gma
      dad.dad = gpa
      uncle.mom = gma
      uncle.dad = gpa
      uncle.kids.add(cousin)
      cousin.dad = uncle
      gma.kids.add(@[dad, uncle])
      mom.kids.add(@[k1, k2])
      dad.kids.add(@[k1, k2])
      ggma.kids.add(gma)
      gma.mom = ggma
      gma.dad = ggpa

      check relatedness(uncle, k1) == 0.25
      check k1.relatedness(cousin) == 0.125


    test "relatedness also":
      var k1 = Sample(family_id:"1", id:"kid1")
      var k2 = Sample(family_id:"1", id:"kid2")
      var mom = Sample(family_id:"1", id:"mom")
      var dad = Sample(family_id:"1", id:"dad")
      var uncle = Sample(family_id: "1", id:"uncle")
      var cousin = Sample(family_id: "1", id:"cousin")
      var gma = Sample(family_id:"1", id:"gma")
      var gpa = Sample(family_id: "1", id:"gpa")
      var ggma = Sample(family_id:"1", id:"ggma")
      var unrel = Sample(family_id:"1", id:"un")
      var extern = Sample(family_id:"xxx", id:"extern")
      k1.mom = mom
      k2.mom = mom
      k1.dad = dad
      k2.dad = dad
      dad.mom = gma
      dad.dad = gpa
      uncle.mom = gma
      uncle.dad = gpa
      uncle.kids.add(cousin)
      cousin.dad = uncle
      gma.kids.add(@[dad, uncle])
      gpa.kids.add(@[dad, uncle])

      mom.kids.add(@[k1, k2])
      dad.kids.add(@[k1, k2])
      ggma.kids.add(gma)
      gma.mom = ggma
      check relatedness(uncle, k1) == 0.25
