#!/bin/bash
#
gfortran -c -g region_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling region_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran region_prb.o -L$HOME/lib/$ARCH -lregion
if [ $? -ne 0 ]; then
  echo "Errors linking and loading region_prb.o"
  exit
fi
rm region_prb.o
#
mv a.out region_prb
./region_prb > region_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running region_prb"
  exit
fi
rm region_prb
#
echo "Test program output written to region_prb_output.txt."
