#!/bin/bash
#
gfortran -c -g asa183_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa183_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa183_prb.o -L$HOME/lib/$ARCH -lasa183
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa183_prb.o"
  exit
fi
rm asa183_prb.o
#
mv a.out asa183_prb
./asa183_prb > asa183_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa183_prb"
  exit
fi
rm asa183_prb
#
echo "Test program output written to asa183_prb_output.txt."
