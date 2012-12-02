#!/bin/bash
#
gfortran -c -g asa147_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa147_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran asa147_prb.o -L$HOME/lib/$ARCH -lasa147
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa147_prb.o"
  exit
fi
rm asa147_prb.o
#
mv a.out asa147_prb
./asa147_prb > asa147_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa147_prb"
  exit
fi
rm asa147_prb
#
echo "Test program output written to asa147_prb_output.txt."
