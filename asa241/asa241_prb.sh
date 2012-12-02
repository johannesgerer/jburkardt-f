#!/bin/bash
#
gfortran -c -g asa241_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa241_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa241_prb.o -L$HOME/lib/$ARCH -lasa241
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa241_prb.o"
  exit
fi
rm asa241_prb.o
#
mv a.out asa241_prb
./asa241_prb > asa241_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa241_prb"
  exit
fi
rm asa241_prb
#
echo "Test program output written to asa241_prb_output.txt."
