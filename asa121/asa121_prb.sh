#!/bin/bash
#
gfortran -c -g asa121_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa121_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa121_prb.o -L$HOME/lib/$ARCH -lasa121
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa121_prb.o"
  exit
fi
rm asa121_prb.o
#
mv a.out asa121_prb
./asa121_prb > asa121_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa121_prb"
  exit
fi
rm asa121_prb
#
echo "Test program output written to asa121_prb_output.txt."
