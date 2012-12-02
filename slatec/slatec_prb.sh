#!/bin/bash
#
gfortran -c -g slatec_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling slatec_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran slatec_prb.o -L$HOME/lib/$ARCH -lslatec
if [ $? -ne 0 ]; then
  echo "Errors linking and loading slatec_prb.o"
  exit
fi
rm slatec_prb.o
#
mv a.out slatec_prb
./slatec_prb > slatec_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running slatec_prb"
  exit
fi
rm slatec_prb
#
echo "Test program output written to slatec_prb_output.txt."
