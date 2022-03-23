{.passC:"-mpopcnt".}

import std/bitops

const seedA:uint64 = 0x3c8bfbb395c60474'u64
const seedC:uint64 = 0x3193c18562a02b4c'u64
const seedG:uint64 = 0x20323ed082572324'u64
const seedT:uint64 = 0x295549f54be24456'u64
const seedN:uint64 = 0x0000000000000000'u64

const offset:uint8 = 0x7

const lookup: array[256, uint64] = [
 seedN, seedT, seedN, seedG, seedA, seedA, seedN, seedC, # 0..7
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 8..15
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 16..23
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 24..31
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 32..39
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 40..47
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 48..55
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 56..63
 seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, # 64..71
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 72..79
 seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, # 80..87
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 88..95
 seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, # 96..103
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 104..111
 seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, # 112..119
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 120..127
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 128..135
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 136..143
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 144..151
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 152..159
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 160..167
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 168..175
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 176..183
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 184..191
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 192..199
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 200..207
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 208..215
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 216..223
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 224..231
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 232..239
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 240..247
 seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, # 248..255
]

    
proc ntf64*(s:string, k:int): uint64 {.inline.} =
  for i in 0..<k:
    result = rotateLeftBits(result, 1'u)
    result = result xor lookup[cast[uint8](s[i])]

proc ntr64(s:string, k:int): uint64 {.inline.} =
  for i in 0..<k:
    result = rotateLeftBits(result, 1'u)
    result = result xor lookup[cast[uint8](s[k - 1 - i]) and offset]


proc forward_nthash*(e:var uint64, remove_base:char, add_base:char, k:int) {.inline.} =
  e = rotateLeftBits(e, 1'u)
  e = e xor rotateLeftBits(lookup[cast[uint8](remove_base)], k)
  e = e xor lookup[cast[uint8](add_base)]

proc reverse_nthash*(e:var uint64, remove_base:char, add_base:char, k:int) {.inline.} =
  e = rotateRightBits(e, 1'u)
  e = e xor rotateRightBits(lookup[cast[uint8](remove_base) and offset], 1)
  e = e xor rotateLeftBits(lookup[cast[uint8](add_base) and offset], k - 1)
 
iterator nthash_forward*(s:string, k:int): uint64 =
  var f = s[0..<k].ntf64(k)
  for i in k..s.high:
    yield f
    f.forward_nthash(s[i - k], s[i], k)
  yield f

iterator nthash*(s:string, k:int): uint64 =
  var f = s[0..<k].ntf64(k)
  var r = s[0..<k].ntr64(k)

  for i in k..s.high:
    yield min(f, r)
    f.forward_nthash(s[i - k], s[i], k)
    r.reverse_nthash(s[i - k], s[i], k)
  yield min(f, r)

#import ./ringbuffer

iterator strobemer_forward*(s:string, k:int, wMin:int, wMax:int): uint64 =
    ## here, k is the `l` from the strobemer paper so the actual bases in the result is 2*k
    ## n is not currently used.
    var h = s[0..<k].ntf64(k)
    var hashes = newSeqOfCap[uint64](s.len)
    hashes.add(h)

    let q:uint64 = (1 shl 23) - 1'u64

    # calculate all hashes and keep in memory. barf
    # TODO: update to use ring-buffer
    for i in 1..<s.len - k + 1:
        h.forward_nthash(s[i - 1], s[k + i - 1], k)
        hashes.add(h)

    # subtracting 2 here so we get at least one shrunk value from inner loop.
    for i in wMax..<s.len + wMax - wMin - 2:
        var h1 = hashes[i - wMax]
        var h2m = uint64.high
        var h2i = i - wMax + wMin
        for j in h2i..<min(i - wMax, hashes.high):
            #if i + j >= len(hashes): break
            # shen  : https://github.com/ksahlin/strobemers/blob/main/randstrobe_implementations/src/index.cpp#L196            
            # sahlin: https://github.com/ksahlin/strobemers/blob/main/randstrobe_implementations/src/index.cpp#L334
            # sahlin popcount: https://github.com/ksahlin/strobemers/blob/main/randstrobe_implementations/src/index.cpp#L335
            #let ht = countSetBits(uint64(h1 xor hashes[j]) and q)
            let ht = countSetBits((h1 xor hashes[j]) and q).uint64
            if ht < h2m:
               h2m = ht
               h2i = j
            
        yield (h1 shr 1) + uint64(float64(hashes[h2i]) / 3)
             
when isMainModule:
  import random
  import strformat
  import sets

  var x = "ACGTGACGTACGT"

  for k in x.nthash_forward(5):
    echo k
  echo ""
  for k in x.nthash(5):
    echo k

  block:
     echo "strobemers"
     # make a random sequence
     randomize(42)
     var s1 = newString(0)
     var chars = {'A', 'C', 'T', 'G'}
     for i in 0..100: s1.add(sample(chars))
     var s2 = deepCopy(s1)
     echo "s1:", s1
     var nmuts = 0
     for i in 0..s2.high:
         if rand(1.0) < 0.05:
            nmuts += 1
            s2[i] = sample(chars - {s2[i]})
     echo "nmuts:", nmuts
     assert s2 != s1
     echo "s2:", s2

     #echo ">", s1
     var t = initHashSet[uint64]()
     var n = 0

     #for k in s1.slide_forward(10):
     for k in strobemer_forward(s1, 3, 3, 7):
         n.inc
         t.incl(k)

     var seen = 0
     #for k in s2.slide_forward(10):
     for k in strobemer_forward(s2, 3, 3, 7):
       if k in t: seen += 1
     echo &"seen:{seen} of:{n} -> {seen.float/n.float*100:.1f}% (higher is better)"
     echo &"unique:{t.len} of:{n} -> {t.len.float/n.float*100:.1f}% (higher is better)"


  if false:
    for str in ["ACGTCGGCGCTTAGCTAGACCACGCTGCACGTC",
                "ACGTCGTCGCTTAGCTAGACCACGCTGCACGTC"]:
     
      echo len(str)
      for k in strobemer_forward(str, 4, 5, 7):
        echo k
      echo "\nfrom:  ", str
