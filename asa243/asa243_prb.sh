#!/bin/bash
#
gfortran -c -g asa243_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa243_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa243_prb.o -L$HOME/lib/$ARCH -lasa243
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa243_prb.o"
  exit
fi
rm asa243_prb.o
#
mv a.out asa243_prb
./asa243_prb > asa243_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa243_prb"
  exit
fi
rm asa243_prb
#
echo "Test program output written to asa243_prb_output.txt."
