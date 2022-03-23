
type RingBuffer*[T] = object
    ## very simple and efficient ring buffer
    data: seq[T]
    head: int

proc newRingBuffer*[T](values: seq[T]): RingBuffer[T] {.noInit.} =
    result = RingBuffer[T](data: values, head: 0)

proc `[]`*[T](b: RingBuffer[T], i: int): T {.inline.} =
    when defined(boundChecks):
        doAssert i + b.head <= 2 * b.data.len - 1
    b.data[(i + b.head) mod b.data.len]

proc add*[T](b: var RingBuffer[T], val: T) {.inline.} =
    b.data[b.head] = val
    b.head.inc
    if b.head == b.data.len: b.head = 0
    
when isMainModule:
    import unittest
    import strformat
    suite "ringbuffer":
      test "add":
          var b = newRingBuffer[uint64](@[1u64, 2, 3, 4, 5])
          check b[0] == 1
          check b[1] == 2
          check b[2] == 3
          check b[4] == 5
          check b[16] == 2

          for i in 11..18:
              echo "i:", $i
              b.add(i.uint64)
              echo b.data
              for i in 0..<b.data.len:
                  stdout.write(&"{i}: {b[i]}, ")
              stdout.writeline("")
