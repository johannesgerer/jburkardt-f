#!/bin/bash
#
gfortran -c -g toeplitz_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toeplitz_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toeplitz_prb.o -L$HOME/lib/$ARCH -ltoeplitz
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toeplitz_prb.o"
  exit
fi
rm toeplitz_prb.o
#
mv a.out toeplitz_prb
./toeplitz_prb > toeplitz_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toeplitz_prb"
  exit
fi
rm toeplitz_prb
#
echo "Test program output written to toeplitz_prb_output.txt."
