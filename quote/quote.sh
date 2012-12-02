#!/bin/bash
#
gfortran -c -g quote.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling quote.f90"
  exit
fi
rm compiler.txt
#
gfortran quote.o
if [ $? -ne 0 ]; then
  echo "Errors while loading quote.o"
  exit
fi
rm quote.o
#
mv a.out ~/bin/$ARCH/quote
#
echo "Executable installed as ~/bin/$ARCH/quote"
