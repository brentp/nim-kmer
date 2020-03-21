## this module provides functions to encode and decode DNA strings (kmers) into uint64's
## to reverse complement the encoded values (in-place), and to add bases to the right or
## left end withtout decoding.
## Any characters other than AaCcTtGg are encoded as A
## The slide iterator yields each encoded (min) kmer along a given input string.

type kmer* = object {.packed.}
  reverse_strand {.bitsize: 1.}: uint64
  mer {.bitsize: 63.}: uint64

# convert a char to its 01234 encoding
const lookup: array[123, uint64] = [1'u64, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1]

proc encode*(k: string, reverse_strand:bool=false): kmer {.inline, noInit.} =
  ## encode a string into a uint64
  when not defined(danger):
    if len(k) > 31:
      quit "[kmer] can't encode kmer longer than 31"
  var tmp:uint64
  for c in k:
    tmp = (tmp or lookup[cast[uint8](c)]) shl 2
  result.reverse_strand = reverse_strand.uint64
  result.mer = tmp shr 2

const str = "CATGN"

proc `~`(k:kmer): kmer {.inline, noInit.} =
  result.mer = not k.mer
  result.reverse_strand = not k.reverse_strand

proc decode*(ek: kmer, k: var string) {.inline.} =
  ## decode a string from a uint64 into k. the length
  ## of k determines how much is decoded
  var
    base {.noInit.}: int
    i = k.len
    L = i
    e:uint64 = ek.mer
  while i > 0:
    base = int((e shr (i * 2  - 2)) and 3'u64)
    k[L-i] = str[base]
    i -= 1

proc decode*(e: kmer, k: var array[6, char], length: int) {.inline.} =
  ## decode a string from a uint64 into k. the length
  ## of k determines how much is decoded
  var
    base {.noInit.}: int
    i = length
    e = e.mer
  while i > 0:
    base = int((e shr (i * 2  - 2)) and 3'u64)
    k[length-i] = str[base]
    i -= 1

proc reverse_complement*(encoded: kmer, L:int|uint64): kmer {.inline, noInit.} =
  ## fast reverse complement of encoded value
  # from Zev: https://www.biostars.org/p/113640/#424280
  var mer = not encoded.mer
  mer = ((mer shr 2'u64 and 0x3333333333333333'u64) or (mer and 0x3333333333333333'u64) shl 2'u64);
  mer = ((mer shr 4'u64 and 0x0F0F0F0F0F0F0F0F'u64) or (mer and 0x0F0F0F0F0F0F0F0F'u64) shl 4'u64);
  mer = ((mer shr 8'u64 and 0x00FF00FF00FF00FF'u64) or (mer and 0x00FF00FF00FF00FF'u64) shl 8'u64);
  mer = ((mer shr 16'u64 and 0x0000FFFF0000FFFF'u64) or (mer and 0x0000FFFF0000FFFF'u64) shl 16'u64);
  mer = ((mer shr 32'u64 and 0x00000000FFFFFFFF'u64) or (mer and 0x00000000FFFFFFFF'u64) shl 32'u64);
  result = kmer(mer:mer shr (2 * (32 - L)), reverse_strand: not encoded.reverse_strand)

proc mincode*(k: string): kmer {.inline, noInit.} =
  ## encode a string into the min of itself and its reverse complement
  let f = k.encode()
  return cast[kmer](min(cast[uint64](f), cast[uint64](f.reverse_complement(k.len))))

proc right*(encoded: kmer, L:uint64): kmer {.inline, noInit.} =
  result.mer = (encoded.mer) and ((1'u64 shl (2'u64*(L-1'u64))) - 1'u64)
  result.reverse_strand = encoded.reverse_strand

proc left*(encoded: kmer, L:uint64): kmer {.inline, noInit.} =
  result.mer = encoded.mer shr 2
  result.reverse_strand = encoded.reverse_strand

proc forward_add*(encoded: var kmer, base: char, L: int) {.inline.} =
  ## drop the first base from an encoded kmer of length L and add a new one.
  ## useful for sliding along a sequence.
  encoded.mer = (encoded.mer shl 2) and uint64((1 shl (2*L)) - 1)
  encoded.mer = (encoded.mer or lookup[cast[int](base)])

proc reverse_add*(rencoded: var kmer, base: char, L: int) {.inline.} =
  ## drop the last base from a rev-comped kmer of length L and add a new
  ## complemented base at the start.
  var mer = (rencoded.mer shr 2)
  var tmp = 3'u64 - lookup[cast[int](base)]
  tmp = tmp shl (L*2-2)
  rencoded.mer = mer or tmp

iterator slide*(s:string, k: int): kmer {.inline.} =
  ## given a string (DNA seq) yield each possible kmer where
  ## the min of the foward and reverse complement is used.
  var f = s[0..<k].encode()
  var r = f.reverse_complement(k)
  for i in k..s.high:
    yield cast[kmer](min(cast[uint64](f), cast[uint64](r)))
    let base:char = s[i]
    f.forward_add(base, k)
    r.reverse_add(base, k)

  yield cast[kmer](min(cast[uint64](f), cast[uint64](r)))

iterator slide_forward*(s:string, k: int): kmer {.inline.} =
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

iterator slide_forward_mask*(s:string, k: int, mask: seq[bool]): kmer {.inline.} =
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

proc reverse_complement_old*(encoded: uint64, L: int): uint64 {.inline.} =
  ## reverse complement an encoded sequence where L is the kmer length.
  var e = encoded
  var base: uint8
  for i in 0..<L:
    base = 3'u8 - uint8(e and 3)
    e = e shr 2
    result = (result or base) shl 2
  result = result shr 2

when isMainModule:
  import random
  import times
  import strutils
  randomize(32)

  var t = cpuTime()

  var k = "CTCCAGCCGGACGCGGCCGGCAGCAGACGCA"
  var s = k
  var e = k.encode()
  e.decode(s)
  if s != k:
    quit "error encoding and decoding"
  var rc = e.reverse_complement(k.len)
  rc.decode(s)
  if s != "TGCGTCTGCTGCCGGCCGCGTCCGGCTGGAG":
       echo "rev comp error", " "
       echo "was:", k
       echo "got:", s, " ", k.len
       echo "exp:TGCGTCTGCTGCCGGCCGCGTCCGGCTGGAG"
       #quit(2)

  e.forward_add('G', k.len)
  e.decode(s)
  if s != "TCCAGCCGGACGCGGCCGGCAGCAGACGCAG":
    echo "error adding base:"
    echo "got ", s
    echo "exp ", "TCCAGCCGGACGCGGCCGGCAGCAGACGCAG"
    #quit(1)

  rc.reverse_add('A', k.len)
  rc.reverse_add('A', k.len)
  rc.decode(s)
  if s != "TTTGCGTCTGCTGCCGGCCGCGTCCGGCTGG":
    echo "error adding reverse base:"
    echo "got ", s
    #quit(1)

  #for i, s in pairs("TTTGCGTCTGCTGCCGGCCGCGTCCGGCTGG"):
  #    echo i, s
  k = k[0..<25]
  s = k

  echo "testing round-trip on random kmers"
  for i in 0..200000:
    shuffle(k)
    e = k.encode()
    e.decode(s)
    var rc = e.mer.reverse_complement_old(s.len)
    var rc2 = e.reverse_complement(s.len)
    if s != k:
       echo "error:", i
       echo "error:", k
       echo "error:", s
       quit(2)
    if rc != rc2.mer:
      echo "error:", i
      echo " k:", k
      echo " s:", s

  echo cpuTime() - t

  echo "testing old reverse_complement"
  let ntimes = 10_000_000
  var n = 0
  var encs = newSeq[uint64](ntimes)
  var kencs = newSeq[kmer](ntimes)
  for i in 0..<ntimes:
    shuffle(k)
    kencs.add(k.encode())
    encs.add(kencs[kencs.high].mer)

  var rco:uint64
  t = cpuTime()
  for e in encs:
    rco = e.reverse_complement_old(s.len)
    if rco >= 65219000000000'u64:
      n.inc
  echo "time:", cpuTime() - t, " n:", n

  echo "testing new reverse_complement"
  t = cpuTime()
  n = 0
  for i, e in kencs:
    rc = e.reverse_complement(s.len)
    if rc.mer >= 65219000000000'u64:
      n.inc
  echo "time:", cpuTime() - t, " n:", n


  var base_str = "CCACGTACTGA"
  var skm = base_str.encode
  var R = skm.right(base_str.len.uint64)
  var L = skm.left(base_str.len.uint64)
  var right_str = newString(base_str.len - 1)
  R.decode(right_str)
  echo "base:", base_str
  echo "right:", right_str
  var left_str = newString(base_str.len - 1)
  L.decode(left_str)
  echo "left:", left_str
  doAssert left_str == base_str[0..<base_str.high]
  doAssert right_str == base_str[1..base_str.high]

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
    (~e).decode(d)
    echo s
    echo d

    echo "before:", e, cast[uint64](e)
    e.reverse_strand = 1
    echo "after:", e, cast[uint64](e)

