## this module provides functions to encode and decode DNA strings (kmers) into uint64's
## to reverse complement the encoded values (in-place), and to add bases to the right or
## left end withtout decoding.
## Any characters other than AaCcTtGg are encoded as A
## The slide iterator yields each encoded (min) kmer along a given input string.

type kmer* = string
## kmer is just a string.

# convert a char to its 01234 encoding
const lookup: array[117, uint64] = [1'u64, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]

proc encode*(k: kmer): uint64 {.inline.} =
  ## encode a string into a uint64
  if len(k) > 31:
    stderr.write_line("[kmer] can't encode kmer longer than 31")
    quit(1)
  for c in k:
    result = (result or lookup[int(c)]) shl 2
  result = result shr 2

const str = "CATGN"

proc decode*(e: uint64, k: var kmer) {.inline.} =
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

proc reverse_complement*(encoded: uint64, L: int): uint64 {.inline.} =
  ## reverse complement an encoded sequence where L is the kmer length.
  var e = encoded
  var base: uint8
  for i in 0..<L:
    base = 3'u8 - uint8(e and 3)
    e = e shr 2
    result = (result or base) shl 2
  result = result shr 2

proc mincode*(k: kmer): uint64 {.inline.} =
  ## encode a string into the min of itself and its reverse complement
  let f = k.encode()
  return min(f, f.reverse_complement(k.len))

proc forward_add*(encoded: var uint64, base: char, L: int) {.inline.} =
  ## drop the first base from an encoded kmer of length L and add a new one.
  ## useful for sliding along a sequence.
  encoded = (encoded shl 2) and uint64((1 shl (2*L)) - 1)
  encoded = (encoded or lookup[int(base)])

proc reverse_add*(rencoded: var uint64, base: char, L: int) {.inline.} =
  ## drop the last base from a rev-comped kmer of length L and add a new
  ## complemented base at the start.
  rencoded = (rencoded shr 2)
  var tmp = 3'u64 - lookup[int(base)]
  tmp = tmp shl (L*2-2)
  rencoded = rencoded or tmp

iterator slide*(s:string, k: int): uint64 {.inline.} =
  ## given a string (DNA seq) yield each possible kmer where
  ## the min of the foward and reverse complement is used.
  var f = s[0..<k].encode()
  var r = f.reverse_complement(k)
  for i in k..s.high:
    yield min(f, r)

    var base:char = s[i]
    f.forward_add(base, k)
    r.reverse_add(base, k)

  yield min(f, r)

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
  echo k.len
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

  for i, s in pairs("TTTGCGTCTGCTGCCGGCCGCGTCCGGCTGG"):
      echo i, s

  echo "testing round-trip on random kmers"
  for i in 0..2000000:
    shuffle(k)
    e = k.encode()
    e.decode(s)
    rc = e.reverse_complement(s.len)
    if s != k:
       echo "error:", i
       echo "error:", k
       echo "error:", s
       quit(2)

  echo cpuTime() - t

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
