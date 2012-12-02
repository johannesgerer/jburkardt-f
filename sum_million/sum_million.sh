#!/bin/bash
#
gfortran -c sum_million.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sum_million.f90"
  exit
fi
rm compiler.txt
#
gfortran sum_million.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sum_million.o"
  exit
fi
rm sum_million.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/sum_million
#
echo "Program installed as ~/bin/$ARCH/sum_million"
