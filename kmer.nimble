# Package

version       = "0.2.7"
author        = "Brent Pedersen"
description   = "dna kmer ops for nim"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.2" #, "nim-lang/c2nim>=0.9.13"
srcDir = "src"


skipDirs = @["tests"]

task test, "run the tests":
  exec "nim c -d:release --lineDir:on --debuginfo -r src/kmer"
