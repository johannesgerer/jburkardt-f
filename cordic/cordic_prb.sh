#!/bin/bash
#
gfortran -c -g cordic_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cordic_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cordic_prb.o -L$HOME/lib/$ARCH -lcordic
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cordic_prb.o"
  exit
fi
rm cordic_prb.o
#
mv a.out cordic_prb
./cordic_prb > cordic_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cordic_prb"
  exit
fi
rm cordic_prb
#
echo "Test program output written to cordic_prb_output.txt."
