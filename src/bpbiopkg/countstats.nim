
type CountStats*[T: SomeOrdinal] = object
  ## store counts of observed values in a fixed-size array
  ## and get median (or any percentile) from those.
  ## we can even keep a moving window of medians as we can drop
  ## values by decrementing the counts.
  counts*: seq[T]
  n*: int
  n_below: int
  current: T
  min_value*: T
  max_size: int32

proc initCountStats*[T](min_value:T=T(0), size:int=128, max_size:int32=262143): CountStats[T] =
  return CountStats[T](counts:newSeq[T](size), min_value:min_value.T, max_size:max_size)

proc adjust*[T](m: var CountStats[T]) {.inline.} =
  var mid = int(0.5 + m.n / 2)
  if mid <= m.n_below:
    while m.current > 0.T:
      m.current -= 1
      m.n_below -= m.counts[m.current.int].int
      if m.counts[m.current.int] != 0: break

  elif mid > m.n_below + m.counts[m.current.int].int:
    while m.current.int < m.counts.len:
      m.n_below += m.counts[m.current.int].int
      m.current += 1
      if m.counts[m.current.int] != 0: break

proc add*[T](m:var CountStats[T], d:T) {.inline.} =
  if d < m.min_value:
    m.counts[0] += 1
  elif d - m.min_value < m.counts.high.T:
    m.counts[d - m.min_value] += 1
  elif m.counts.len < m.max_size:
    var o = m.counts.len
    m.counts.setLen(min(m.max_size, (1.5*m.counts.len.float).int))
    zeroMem(m.counts[o].addr, (m.counts.len - o) * sizeof(T))
    m.counts[min(m.counts.high.T, d - m.min_value)] += 1
  else:
    m.counts[m.counts.high.T] += 1
  if d.int - m.min_value.int < m.current.int:
    m.n_below += 1

  m.n += 1
  m.adjust()

proc drop*[T](m:var CountStats[T], d:T) {.inline.} =
  if d < m.min_value:
    m.counts[0] -= 1
  elif d - m.min_value < m.counts.high.T:
    m.counts[d - m.min_value] -= 1
  else:
    m.counts[m.counts.high] -= 1
  if d.int - m.min_value.int < m.current.int:
    m.n_below -= 1
  m.n -= 1
  m.adjust()


proc median*[T](m:var CountStats[T]): T {.inline.} =
  return m.current + m.min_value

proc mean*[T](m:var CountStats[T]): float64 {.inline.} =
  for i, cnt in m.counts:
    result += (i * cnt)
  result = (0.5 + result.float64 / m.n.float64) + m.min_value.float64

proc clear*[T](m:var CountStats[T]) =
  zeroMem(m.counts[0].addr, sizeof(m.counts[0]) * m.counts.len)
  m.current = T(0)
  m.n = 0
  m.n_below = 0


proc percentile*[T](m:CountStats[T], pct:float64): int {.inline.} =
  ## return the value at the requested percentile.
  var cum = 0
  var p = pct
  if p > 1: p /= 100.0
  if p > 1:
      quit "can't get value outside of distribution (> 1)"

  var stop_n = (m.n.float64 * p).int
  for i, cnt in m.counts:
    cum += cnt
    if cum >= stop_n:
      return i + m.min_value
  return -1

proc `$`*[T](m:var CountStats[T]): string =
  return &"CountStats(n:{m.n}, median:{m.median} vals: {m.counts[0..100]})"

when isMainModule:
    var c = CountStats[int](counts:newSeq[int](256))
    c.add(1)
    assert c.median == 1
    c.add(3)
    c.add(2)
    assert c.median == 2

    import unittest

    suite "count test suite":
      test "test median":

          var c = CountStats[int](counts:newSeq[int](256))
          c.add(1)
          check c.median == 1

          c.add(3)
          c.add(2)
          check c.median == 2
          check c.n == 3

          c.drop(2)
          c.drop(3)
          check c.median == 1
          check c.n == 1

          for i in 0..<100:
            c.add(100)
          check c.n == 101
          check c.median == 100

      test "init":
        var c = initCountStats[uint8](min_value=uint8(32))
        c.add(33)
        check c.median == 33
        c.add(35)
        c.add(34)
        check c.median == 34
        c.drop(35)
        c.drop(34)
        check c.median == 33


      test "example":
        var c = initCountStats[uint16](min_value=uint16(328))
        c.add(400'u16)
        c.add(401'u16)
        c.add(403'u16)
        check c.median == 401'u16
        for i in 0..100:
          c.add(401'u16)
          check c.median == 401'u16
          check c.n == 4 + i

        c.drop(400)
        check c.median == 401
        c.add(401)
        check c.median == 401

        for j in 0..100:
          c.drop(401'u16)
          check c.median == 401'u16

        c.drop(401)
        c.drop(401)
        echo c.median
