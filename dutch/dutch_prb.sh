#!/bin/bash
#
gfortran -c -g dutch_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dutch_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran dutch_prb.o -L$HOME/lib/$ARCH -ldutch
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dutch_prb.o"
  exit
fi
rm dutch_prb.o
#
mv a.out dutch_prb
./dutch_prb > dutch_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dutch_prb"
  exit
fi
rm dutch_prb
#
echo "Test program output written to dutch_prb_output.txt."
