#!/bin/bash
#
gfortran -c -g latin_center_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_center_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran latin_center_prb.o -L$HOME/lib/$ARCH -llatin_center
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_center_prb.o"
  exit
fi
rm latin_center_prb.o
#
mv a.out latin_center_prb
./latin_center_prb > latin_center_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latin_center_prb"
  exit
fi
rm latin_center_prb
#
echo "Test program output written to latin_center_prb_output.txt."
