#!/bin/bash
#
gfortran -c -g colors_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling colors_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran colors_prb.o -L$HOME/lib/$ARCH -lcolors
if [ $? -ne 0 ]; then
  echo "Errors linking and loading colors_prb.o"
  exit
fi
rm colors_prb.o
#
mv a.out colors_prb
./colors_prb > colors_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running colors_prb"
  exit
fi
rm colors_prb
#
echo "Test program output written to colors_prb_output.txt."
