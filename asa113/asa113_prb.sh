#!/bin/bash
#
gfortran -c -g asa113_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa113_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa113_prb.o -L$HOME/lib/$ARCH -lasa113
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa113_prb.o"
  exit
fi
rm asa113_prb.o
#
mv a.out asa113_prb
./asa113_prb > asa113_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa113_prb"
  exit
fi
rm asa113_prb
#
echo "Test program output written to asa113_prb_output.txt."
