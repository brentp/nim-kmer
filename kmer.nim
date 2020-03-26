## this module provides functions to encode and decode DNA strings (kmers) into uint64's
## to reverse complement the encoded values (in-place), and to add bases to the right or
## left end withtout decoding.
## Any characters other than AaCcTtGg are encoded as A
## The slide iterator yields each encoded (min) kmer along a given input string.

# convert a char to its 01234 encoding
const lookup: array[123, uint64] = [1'u64, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1]

proc encode*(k: string): uint64 {.inline.} =
  ## encode a string into a uint64
  when not defined(danger):
    if len(k) > 31:
      quit "[kmer] can't encode kmer longer than 31"
  for c in k:
    result = (result or lookup[cast[uint8](c)]) shl 2
  result = result shr 2

const str = "CATGN"

proc decode*(e: uint64, k: var string) {.inline.} =
  ## decode a string from a uint64 into k. the length
  ## of k determines how much is decoded
  var
    base: int
    i = k.len
    L = i
  while i > 0:
    base = int((e shr (i * 2  - 2)) and 3'u64)
    k[L-i] = str[base]
    i -= 1

proc reverse_complement*(encoded: uint64, L:int|uint64): uint64 {.inline.} =
  ## fast reverse complement of encoded value
  # from Zev: https://www.biostars.org/p/113640/#424280
  result = not encoded
  result = ((result shr 2'u64 and 0x3333333333333333'u64) or (result and 0x3333333333333333'u64) shl 2'u64);
  result = ((result shr 4'u64 and 0x0F0F0F0F0F0F0F0F'u64) or (result and 0x0F0F0F0F0F0F0F0F'u64) shl 4'u64);
  result = ((result shr 8'u64 and 0x00FF00FF00FF00FF'u64) or (result and 0x00FF00FF00FF00FF'u64) shl 8'u64);
  result = ((result shr 16'u64 and 0x0000FFFF0000FFFF'u64) or (result and 0x0000FFFF0000FFFF'u64) shl 16'u64);
  result = ((result shr 32'u64 and 0x00000000FFFFFFFF'u64) or (result and 0x00000000FFFFFFFF'u64) shl 32'u64);
  return (result shr (2 * (32 - L)));

proc forward_add*(encoded: var uint64, base: char, L: int) {.inline.} =
  ## drop the first base from an encoded kmer of length L and add a new one.
  ## useful for sliding along a sequence.
  encoded = (encoded shl 2) and ((1'u64 shl (2*L)) - 1)
  encoded = (encoded or lookup[cast[int](base)])

proc reverse_add*(rencoded: var uint64, base: char, L: int) {.inline.} =
  ## drop the last base from a rev-comped kmer of length L and add a new
  ## complemented base at the start.
  rencoded = (rencoded shr 2)
  var tmp = 3'u64 - lookup[cast[int](base)]
  tmp = tmp shl (L*2-2)
  rencoded = rencoded or tmp

type stranded* = tuple[enc:uint64, min_complement:uint8]

proc mincode*(k: string): stranded {.inline.} =
  ## encode a string into the min of itself and its reverse complement
  let f = k.encode()
  let r = f.reverse_complement(k.len)
  return (min(f, r), cast[uint8](r < f))

iterator slide*(s:string, k: int): stranded {.inline.} =
  ## given a string (DNA seq) yield each possible kmer where
  ## the min of the foward and reverse complement is used.
  var f = s[0..<k].encode()
  var r = f.reverse_complement(k)
  var base:char
  for i in k..s.high:
    yield (min(f, r), cast[uint8](r < f))
    base = s[i]
    f.forward_add(base, k)
    r.reverse_add(base, k)
  yield (min(f, r), cast[uint8](r < f))


iterator slide_forward*(s:string, k: int): uint64 {.inline.} =
  ## given a string (DNA seq) yield each possible kmer where
  ## the min of the foward and reverse complement is used.
  var f = s[0..<k].encode()
  for i in k..s.high:
    yield f
    let base:char = s[i]
    f.forward_add(base, k)
  yield f


proc sum(m: seq[bool]): int {.inline.} =
  for v in m: result += v.int

iterator slide_forward_mask*(s:string, k: int, mask: seq[bool]): uint64 {.inline.} =
  ## take a boolean mask choosing which bases to extract. this is less
  ## efficient than slide_forward, bt allows sparse kmers
  assert k == sum(mask)
  var sequence = newStringOfCap(k)
  for i, v in mask:
    if v: sequence.add(s[i])

  var f = sequence.encode()
  yield f
  for j in 1..s.high - mask.high:
    for i, v in mask:
      if v: f.forward_add(s[i + j], k)
    yield f


iterator dists*(s: string, k:int): auto {.inline.} =
  ## yield each (min) encoded k-mer and its distance from the closest end of the read.
  var i = 0
  var m = s.len - k
  for x in s.slide(k):
    yield (min(i, m - i), x)
    inc i


when isMainModule:
  import random
  import times
  import strutils

  var t = cpuTime()

  var k = "CTCCAGCCGGACGCGGCCGGCAGCAGACGCA"
  var s = k
  var e = k.encode()
  e.decode(s)
  var rc = e.reverse_complement(k.len)
  rc.decode(s)
  if s != "TGCGTCTGCTGCCGGCCGCGTCCGGCTGGAG":
       echo "rev comp error"
       quit(2)

  e.forward_add('G', k.len)
  e.decode(s)
  if s != "TCCAGCCGGACGCGGCCGGCAGCAGACGCAG":
    echo "error adding base:"
    echo "got ", s
    quit(1)

  rc.reverse_add('A', k.len)
  rc.reverse_add('A', k.len)
  rc.decode(s)
  if s != "TTTGCGTCTGCTGCCGGCCGCGTCCGGCTGG":
    echo "error adding reverse base:"
    echo "got ", s
    quit(1)

  #for i, s in pairs("TTTGCGTCTGCTGCCGGCCGCGTCCGGCTGG"):
  #    echo i, s
  k = k[0..<25]
  s = k
  proc reverse_complement_old*(encoded: uint64, L: int): uint64 {.inline.} =
    ## reverse complement an encoded sequence where L is the kmer length.
    var e = encoded
    var base: uint8
    for i in 0..<L:
      base = 3'u8 - uint8(e and 3)
      e = e shr 2
      result = (result or base) shl 2
    result = result shr 2

  echo "testing round-trip on random kmers"
  for i in 0..200000:
    shuffle(k)
    e = k.encode()
    e.decode(s)
    rc = e.reverse_complement_old(s.len)
    var rc2 = e.reverse_complement(s.len)
    if s != k:
       echo "error:", i
       echo "error:", k
       echo "error:", s
       quit(2)
    if rc != rc2:
      echo "error:", i, " k:", k, " s:", s

  echo cpuTime() - t

  echo "testing old reverse_complement"
  let ntimes = 10_000_000
  var n = 0
  var encs = newSeq[uint64](ntimes)
  for i in 0..<ntimes:
    shuffle(k)
    encs.add(k.encode())

  var old = newSeqOfCap[uint64](ntimes)
  t = cpuTime()
  for e in encs:
    rc = e.reverse_complement_old(s.len)
    if rc >= 65219000000000'u64:
      n.inc
  # old.add(rc)
  echo "time:", cpuTime() - t, " n:", n

  echo "testing new reverse_complement"
  t = cpuTime()
  n = 0
  for i, e in encs:
    rc = e.reverse_complement(s.len)
    if rc >= 65219000000000'u64:
      n.inc
  #doAssert old[i] == rc
  echo "time:", cpuTime() - t, " n:", n


  block:
    var s = "ACTGGC"
  #[

  s = "TTTGCGTCTGCTGCCGGCCGCGTCCGGCTGG"
  var space = ""
  var kx1 = new_string(10)
  var kx2 = new_string(10)
  var i = 0
  for d, u in s.dists(10):
    echo s
    u.decode(kx1)
    u.reverse_complement(10).decode(kx2)
    assert kx1 == s[i..<(i+10)] or kx2 == s[i..<(i+10)]
    if kx1 == s[i..<(i+10)]:
      echo space & kx1, " ", d
    else:
      echo space & kx2, " ", d
    space &= " "
    echo ""
    i += 1
  ]#


  block:

    var s = "ACACACACACT"
    var k = 3
    var mask = @[true, false, true, false, true]
    echo "slide mask:", s, " mask:", mask
    for v in s.slide_forward_mask(k, mask):
      var n = newString(k)
      v.decode(n)
      echo n

    mask = @[true, true, false, false, false, true]
    echo "mask:", mask
    for v in s.slide_forward_mask(k, mask):
      var n = newString(k)
      v.decode(n)
      echo n



  block:
    var s = "ACACACACACT"
    var e = s.encode
    var d = newString(s.len)
    (not e).decode(d)
    echo s
    echo d


  block:
    var s = "ACTGACGGACCCGAGGGCACCCGAGGCCTTTTTTTTGCGGGAGGAGGAGACTGACTGCGGGAGGAGGAGACTGACTGCGGGAGGAGGAGACTGACTGCGGGAGGAGGAG"
    var t = cpuTime()

    var lastS = 0'u64
    for i in 0..600_000:
      var S = 0'u64
      for k in s.slide(25):
        S += k.enc mod 15
      if i > 0:
        doAssert S == lastS
      else:
        lastS = S
        echo lastS
        var lastS = S

    echo "slide time:",  cpuTime() - t
