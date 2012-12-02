#!/bin/bash
#
gfortran -c -g spaeth_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spaeth_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran spaeth_prb.o -L$HOME/lib/$ARCH -lspaeth
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spaeth_prb.o"
  exit
fi
rm spaeth_prb.o
#
mv a.out spaeth_prb
#
#  Copy the data files.
#
cp ../../datasets/spaeth/spaeth_*.txt .
#
./spaeth_prb > spaeth_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running spaeth_prb"
  exit
fi
rm spaeth_prb
#
#  Discard the data sets.
#
rm spaeth_*.txt
#
echo "Test program output written to spaeth_prb_output.txt."
