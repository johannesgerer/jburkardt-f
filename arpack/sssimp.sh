#!/bin/bash
#
gfortran -c -g sssimp.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sssimp.f90"
  exit
fi
rm compiler.txt
#
gfortran sssimp.o -L$HOME/lib/$ARCH -larpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sssimp.o"
  exit
fi
rm sssimp.o
#
mv a.out sssimp
./sssimp > sssimp_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sssimp"
  exit
fi
rm sssimp
#
echo "Test program output written to sssimp_output.txt."
