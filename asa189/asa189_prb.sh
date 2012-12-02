#!/bin/bash
#
gfortran -c -g asa189_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa189_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa189_prb.o -L$HOME/lib/$ARCH -lasa189
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa189_prb.o"
  exit
fi
rm asa189_prb.o
#
mv a.out asa189_prb
./asa189_prb > asa189_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa189_prb"
  exit
fi
rm asa189_prb
#
echo "Test program output written to asa189_prb_output.txt."
