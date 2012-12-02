#!/bin/bash
#
gfortran -c -g filum_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filum_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran filum_prb.o -L$HOME/lib/$ARCH -lfilum
if [ $? -ne 0 ]; then
  echo "Errors linking and loading filum_prb.o"
  exit
fi
rm filum_prb.o
#
mv a.out filum_prb
./filum_prb > filum_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running filum_prb"
  exit
fi
rm filum_prb
#
echo "Test program output written to filum_prb_output.txt."
