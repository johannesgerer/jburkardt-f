#!/bin/bash
#
gfortran -c -g faure_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling faure_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran faure_prb.o -L$HOME/lib/$ARCH -lfaure
if [ $? -ne 0 ]; then
  echo "Errors linking and loading faure_prb.o"
  exit
fi
rm faure_prb.o
#
mv a.out faure_prb
./faure_prb > faure_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running faure_prb"
  exit
fi
rm faure_prb
#
echo "Test program output written to faure_prb_output.txt."
