#!/bin/bash
#
gfortran -c -g bvls_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bvls_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran bvls_prb.o $HOME/lib/$ARCH/bvls.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bvls_prb.o + bvls.o"
  exit
fi
rm bvls_prb.o
#
mv a.out bvls_prb
./bvls_prb > bvls_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bvls_prb"
  exit
fi
rm bvls_prb
#
echo "Test program output written to bvls_prb_output.txt."
