#!/bin/bash
#
gfortran -c -g van_der_corput_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling van_der_corput_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran van_der_corput_prb.o -L$HOME/lib/$ARCH -lvan_der_corput
if [ $? -ne 0 ]; then
  echo "Errors linking and loading van_der_corput_prb.o"
  exit
fi
rm van_der_corput_prb.o
#
mv a.out van_der_corput_prb
./van_der_corput_prb > van_der_corput_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running van_der_corput_prb"
  exit
fi
rm van_der_corput_prb
#
echo "Test program output written to van_der_corput_prb_output.txt."
