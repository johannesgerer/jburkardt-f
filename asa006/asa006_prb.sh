#!/bin/bash
#
gfortran -c -g asa006_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa006_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa006_prb.o -L$HOME/lib/$ARCH -lasa006
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa006_prb.o"
  exit
fi
rm asa006_prb.o
#
mv a.out asa006_prb
./asa006_prb > asa006_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa006_prb"
  exit
fi
rm asa006_prb
#
echo "Test program output written to asa006_prb_output.txt."
