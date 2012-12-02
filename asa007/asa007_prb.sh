#!/bin/bash
#
gfortran -c -g asa007_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa007_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa007_prb.o -L$HOME/lib/$ARCH -lasa007
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa007_prb.o"
  exit
fi
rm asa007_prb.o
#
mv a.out asa007_prb
./asa007_prb > asa007_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa007_prb"
  exit
fi
rm asa007_prb
#
echo "Test program output written to asa007_prb_output.txt."
