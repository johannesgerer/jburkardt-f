#!/bin/bash
#
gfortran -c -g cell_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cell_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cell_prb.o -L$HOME/lib/$ARCH -lcell
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cell_prb.o"
  exit
fi
rm cell_prb.o
#
mv a.out cell_prb
./cell_prb > cell_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cell_prb"
  exit
fi
rm cell_prb
#
echo "Test program output written to cell_prb_output.txt."
