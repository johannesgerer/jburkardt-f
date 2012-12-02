#!/bin/bash
#
gfortran -c -g asa091_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa091_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa091_prb.o -L$HOME/lib/$ARCH -lasa091
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa091_prb.o"
  exit
fi
rm asa091_prb.o
#
mv a.out asa091_prb
./asa091_prb > asa091_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa091_prb"
  exit
fi
rm asa091_prb
#
echo "Test program output written to asa091_prb_output.txt."
