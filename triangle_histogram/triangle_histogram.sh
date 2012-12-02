#!/bin/bash
#
gfortran -c -g triangle_histogram.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_histogram.f90"
  exit
fi
rm compiler.txt
#
gfortran triangle_histogram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_histogram.o"
  exit
fi
rm triangle_histogram.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangle_histogram
#
echo "Executable installed as ~/bin/$ARCH/triangle_histogram"
