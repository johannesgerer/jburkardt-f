#!/bin/bash
#
gfortran -c -g table_histogram.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_histogram.f90"
  exit
fi
rm compiler.txt
#
gfortran table_histogram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_histogram.o"
  exit
fi
rm table_histogram.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_histogram
#
echo "Executable installed as ~/bin/$ARCH/table_histogram"
