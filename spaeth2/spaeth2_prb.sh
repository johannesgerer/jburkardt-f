#!/bin/bash
#
gfortran -c -g spaeth2_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spaeth2_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran spaeth2_prb.o -L$HOME/lib/$ARCH -lspaeth2
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spaeth2_prb.o"
  exit
fi
rm spaeth2_prb.o
#
cp ../../datasets/spaeth2/spaeth2_*.txt .
#
mv a.out spaeth2_prb
./spaeth2_prb > spaeth2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running spaeth2_prb"
  exit
fi
rm spaeth2_prb
#
rm spaeth2_0*.txt
rm spaeth2_1*.txt
#
echo "Test program output written to spaeth2_prb_output.txt."
