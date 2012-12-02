#!/bin/bash
#
gfortran -c -g errors_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling errors_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran errors_prb.o -L$HOME/lib/$ARCH -lerrors
if [ $? -ne 0 ]; then
  echo "Errors linking and loading errors_prb.o"
  exit
fi
rm errors_prb.o
#
mv a.out errors_prb
./errors_prb > errors_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running errors_prb"
  exit
fi
rm errors_prb
#
echo "Test program output written to errors_prb_output.txt."
