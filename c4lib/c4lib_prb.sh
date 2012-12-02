#!/bin/bash
#
gfortran -c -g c4lib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c4lib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran c4lib_prb.o -L$HOME/lib/$ARCH -lc4lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c4lib_prb.o"
  exit
fi
rm c4lib_prb.o
#
mv a.out c4lib_prb
./c4lib_prb > c4lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c4lib_prb"
  exit
fi
rm c4lib_prb
#
echo "Test program output written to c4lib_prb_output.txt."
