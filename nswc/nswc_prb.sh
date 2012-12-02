#!/bin/bash
#
gfortran -c -g nswc_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nswc_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran nswc_prb.o -L$HOME/lib/$ARCH -lnswc
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nswc_prb.o"
  exit
fi
rm nswc_prb.o
#
mv a.out nswc_prb
./nswc_prb > nswc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nswc_prb"
  exit
fi
rm nswc_prb
#
echo "Program output written to nswc_prb_output.txt"
