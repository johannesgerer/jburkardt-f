#!/bin/bash
#
gfortran -c -g triangulation_histogram.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_histogram.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_histogram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_histogram.o"
  exit
fi
rm triangulation_histogram.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_histogram
#
echo "Executable installed as ~/bin/$ARCH/triangulation_histogram"
