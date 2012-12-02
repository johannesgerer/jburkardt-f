#!/bin/bash
#
gfortran -c -g asa226_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa226_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa226_prb.o -L$HOME/lib/$ARCH -lasa226
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa226_prb.o"
  exit
fi
rm asa226_prb.o
#
mv a.out asa226_prb
./asa226_prb > asa226_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa226_prb"
  exit
fi
rm asa226_prb
#
echo "Test program output written to asa226_prb_output.txt."
