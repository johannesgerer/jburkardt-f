#!/bin/bash
#
gfortran -c -g wavelet_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wavelet_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran wavelet_prb.o -L$HOME/lib/$ARCH -lwavelet
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wavelet_prb.o"
  exit
fi
rm wavelet_prb.o
#
mv a.out wavelet_prb
./wavelet_prb > wavelet_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wavelet_prb"
  exit
fi
rm wavelet_prb
#
echo "Test program output written to wavelet_prb_output.txt."
