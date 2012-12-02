#!/bin/bash
#
gfortran -c -g asa005_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa005_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa005_prb.o -L$HOME/lib/$ARCH -lasa005
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa005_prb.o"
  exit
fi
rm asa005_prb.o
#
mv a.out asa005_prb
./asa005_prb > asa005_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa005_prb"
  exit
fi
rm asa005_prb
#
echo "Test program output written to asa005_prb_output.txt."
