#!/bin/bash
#
gfortran -c file_row_reverse.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while compiling file_row_reverse.f90"
  exit
fi
rm compiler.txt
#
gfortran file_row_reverse.o
if [ $? -ne 0 ]; then
  echo "Errors occurred while linking file_row_reverse.o"
  exit
fi
rm file_row_reverse.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/file_row_reverse
#
echo "Program installed as ~/bin/$ARCH/file_row_reverse"
