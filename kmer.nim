## this module provides functions to encode and decode DNA strings (kmers) into uint64's
## to reverse complement the encoded values (in-place), and to add bases to the right or
## left end withtout decoding.
## Any characters other than AaCcTtGg are encoded as A
## The slide iterator yields each encoded (min) kmer along a given input string.

# convert a char to its 01234 encoding
#const lookup: array[123, uint64] = [1'u64, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#       1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#       1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 3,
#       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1]


const lookup: array[256, uint64] = [
  0'u64, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
]

const reverse_lookup* = "ACGTN"

type kmer* = array[2, uint64]

proc `<`*(f, r:kmer): bool {.inline.} =
  return f[0] < r[0] or (f[0] == r[0] and f[1] < r[1])

proc min*(f, r:kmer): kmer {.inline.} =
  if f[0] < r[0]: return f
  if r[0] < f[0]: return r
  if f[1] <= r[1]: return f
  return r

template mask(k:int): uint64 =
  (1'u64 shl k) - 1

proc forward_add*(encoded: var kmer, base: char, L: int) {.inline.} =
  ## drop the first base from an encoded kmer of length L and add a new one.
  ## useful for sliding along a sequence.
  let mask = L.mask
  let c = lookup[cast[uint8](base)]
  encoded[0] = ((encoded[0] shl 1) or (c and 1)) and mask
  encoded[1] = ((encoded[1] shl 1) or (c shr 1)) and mask

proc encode*(k: string): kmer {.inline.} =
  ## encode a string into a uint64
  let mask = k.len.mask
  when not defined(danger):
    if len(k) > 63:
      quit "[kmer] can't encode kmer longer than 31"
  for base in k:
    let c = lookup[cast[uint8](base)]
    result[0] = ((result[0] shl 1) or (c and 1)) and mask
    result[1] = ((result[1] shl 1) or (c shr 1)) and mask

proc reverse_add*(rencoded: var kmer, base: char, L: int) {.inline.} =
  ## drop the last base from a rev-comped kmer of length L and add a new
  ## complemented base at the start.
  let c = lookup[cast[uint8](base)]
  rencoded[0] = (rencoded[0] shr 1) or (1 xor (c and 1)) shl (L - 1)
  rencoded[1] = (rencoded[1] shr 1) or (1 xor c shr 1) shl (L - 1)

proc rencode*(k: string): kmer {.inline.} =
  ## reverse encode a string into a uint64
  when not defined(danger):
    if len(k) > 63:
      quit "[kmer] can't encode kmer longer than 31"
  let km1 = uint64(k.len - 1)
  for base in k:
    let c = lookup[cast[uint8](base)]
    result[0] = (result[0] shl 1) or (1 xor (c and 1)) shl km1
    result[1] = (result[1] shl 1) or (1 xor c shr 1) shl km1

type stranded* = tuple[enc:kmer, min_complement:uint8]

proc mincode*(k: string): stranded {.inline.} =
  ## encode a string into the min of itself and its reverse complement
  let f = k.encode()
  let r = k.rencode()
  if f < r:
    return (f, 0'u8)
  return (r, 1'u8)

proc decode*(e: kmer, k: var string) {.inline.} =
  ## decode a string from a uint64 into k. the length
  ## of k determines how much is decoded
  for l in 0..<k.len:
    k[k.len - 1 - l] = reverse_lookup[((e[1] shr l and 1) shl 1) or (e[0] shr l and 1)]

iterator slide*(s:string, k: int): stranded {.inline.} =
  ## given a string (DNA seq) yield each possible kmer where
  ## the min of the foward and reverse complement is used.
  var f = s[0..<k].encode()
  var r = s[0..<k].rencode()
  var base:char
  for i in k..s.high:
    yield (min(f, r), cast[uint8](r < f))
    base = s[i]
    f.forward_add(base, k)
    r.reverse_add(base, k)
  yield (min(f, r), cast[uint8](r < f))

proc sum(m: seq[bool]): int {.inline.} =
  for v in m: result += v.int

iterator slide_forward*(s:string, k: int): kmer {.inline.} =
  ## given a string (DNA seq) yield each possible kmer where
  ## the min of the foward and reverse complement is used.
  var f = s[0..<k].encode()
  for i in k..s.high:
    yield f
    let base:char = s[i]
    f.forward_add(base, k)
  yield f

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


when isMainModule:
  import random
  import times
  import strutils

  var t = cpuTime()

  var k = "CTCCAGCCGGACGCGGCCGGCAGCAGACGCACTCCAGCCGGACGCGGCCGGCAGCAGACGCA"
  var s = k
  var e = k.encode()
  e.decode(s)

  echo k == s

  echo k, " -> ", s

  echo "add 'T'"
  e.forward_add('T', k.len)
  e.decode(s)
  echo s

  echo "add reverse 'T'"
  e = k.encode()
  e.reverse_add('T', k.len)
  e.decode(s)
  echo s

  s.setLen(13)
  echo "s:", k
  for km in k.slide(13):
    km.enc.decode(s)
    echo s, " ", km.min_complement

