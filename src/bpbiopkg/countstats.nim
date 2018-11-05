
type CountStats*[T: SomeOrdinal] = object
  ## store counts of observed values in a fixed-size array
  ## and get median (or any percentile) from those.
  ## we can even keep a moving window of medians as we can drop
  ## values by decrementing the counts.
  counts*: seq[T]
  i_median: T
  n*: int
  i_ok: bool

proc len*[T](m:CountStats[T]): int =
  return m.n

proc add*[T](m:var CountStats[T], d:T) {.inline.} =
  if d < m.counts.high.T:
    m.counts[d] += 1
  else:
    m.counts[m.counts.high.T] += 1
  m.n += 1
  m.i_ok = false

proc drop*[T](m:var CountStats[T], d:T) {.inline.} =
  if d < m.counts.high:
    m.counts[d] -= 1
  else:
    m.counts[m.counts.high] -= 1
  m.n -= 1
  m.i_ok = false

proc median*[T](m:var CountStats[T], skip_zeros:bool=false): T {.inline.} =
  if m.i_ok:
    return m.i_median

  var cum = 0
  var stop_n = (0.5 + m.n.float64 * 0.5).int
  if skip_zeros:
    stop_n = (0.5 + float64(m.n - m.counts[0]) * 0.5).int

  for i, cnt in m.counts:
    if skip_zeros and i == 0: continue
    cum += cnt
    if cum >= stop_n:
      m.i_median = i
      m.i_ok = true
      break
  return m.i_median

proc mean*[T](m:var CountStats[T]): float64 {.inline.} =
  for i, cnt in m.counts:
    result += (i * cnt)
  result = (0.5 + result.float64 / m.n.float64)

proc clear*[T](m:var CountStats[T]) =
  zeroMem(m.counts[0].addr, sizeof(m.counts[0]) * m.counts.len)
  m.i_ok = false
  m.i_median = 0
  m.n = 0

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
      return i
  return -1

proc `$`*[T](m:var CountStats[T]): string =
  return &"CountStats(n:{m.n}, median:{m.median} vals: {m.counts[0..100]})"

when isMainModule:

    import unittest

    suite "count test suite":
      test "test median":

          var c = CountStats[int](counts:newSeq[int](256))
          c.add(1)
          check c.median == 1

          c.add(3)
          c.add(2)
          check c.median == 2
          check c.len == 3

          c.drop(2)
          c.drop(3)
          check c.median == 1
          check c.len == 1

          for i in 0..<100:
            c.add(100)
          check c.len == 101
          check c.median == 100

