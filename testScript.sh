#!/bin/bash
echo "Two-dimensional checks"

echo "Checking fast-marching implementation"
./redistMain 2  32  32 1
./redistMain 2  64  64 1
./redistMain 2 128 128 1
./redistMain 2 256 256 1

echo "Checking quadratic convergence"
./redistMain 2  32  32 2
./redistMain 2  64  64 2
./redistMain 2 128 128 2
./redistMain 2 256 256 2

echo "Checking cubic convergence"
./redistMain 2  32  32 3
./redistMain 2  64  64 3
./redistMain 2 128 128 3
./redistMain 2 256 256 3

echo "Three-dimensional checks"

echo "Checking fast-marching implementation"
./redistMain 3  16  16 1
./redistMain 3  32  32 1
./redistMain 3  64  64 1

echo "Checking quadratic convergence"
./redistMain 3  16  16 2
./redistMain 3  32  32 2
./redistMain 3  64  64 2

#echo "Checking cubic convergence"   <-- not implemented!
#./redistMain 3  16  16 3
#./redistMain 3  32  32 3
#./redistMain 3  64  64 3
