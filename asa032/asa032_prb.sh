#!/bin/bash
#
gfortran -c -g asa032_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa032_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa032_prb.o -L$HOME/lib/$ARCH -lasa032
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa032_prb.o"
  exit
fi
rm asa032_prb.o
#
mv a.out asa032_prb
./asa032_prb > asa032_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa032_prb"
  exit
fi
rm asa032_prb
#
echo "Test program output written to asa032_prb_output.txt."
