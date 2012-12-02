#!/bin/bash
#
gfortran -c -g asa103_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa103_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa103_prb.o -L$HOME/lib/$ARCH -lasa103
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa103_prb.o"
  exit
fi
rm asa103_prb.o
#
mv a.out asa103_prb
./asa103_prb > asa103_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa103_prb"
  exit
fi
rm asa103_prb
#
echo "Test program output written to asa103_prb_output.txt."
