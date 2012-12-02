#!/bin/bash
#
gfortran -c -g asa076_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa076_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa076_prb.o -L$HOME/lib/$ARCH -lasa076
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa076_prb.o"
  exit
fi
rm asa076_prb.o
#
mv a.out asa076_prb
./asa076_prb > asa076_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa076_prb"
  exit
fi
rm asa076_prb
#
echo "Test program output written to asa076_prb_output.txt."
