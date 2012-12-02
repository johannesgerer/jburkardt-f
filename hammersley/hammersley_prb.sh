#!/bin/bash
#
gfortran -c -g hammersley_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hammersley_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran hammersley_prb.o -L$HOME/lib/$ARCH -lhammersley
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hammersley_prb.o"
  exit
fi
rm hammersley_prb.o
#
mv a.out hammersley_prb
./hammersley_prb > hammersley_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hammersley_prb"
  exit
fi
rm hammersley_prb
#
echo "Test program output written to hammersley_prb_output.txt."
