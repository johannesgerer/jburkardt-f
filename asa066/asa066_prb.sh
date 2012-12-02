#!/bin/bash
#
gfortran -c -g asa066_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa066_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa066_prb.o -L$HOME/lib/$ARCH -lasa066
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa066_prb.o"
  exit
fi
rm asa066_prb.o
#
mv a.out asa066_prb
./asa066_prb > asa066_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa066_prb"
  exit
fi
rm asa066_prb
#
echo "Test program output written to asa066_prb_output.txt."
