#!/bin/bash
#
gfortran -c file_transpose.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while compiling file_transpose.f90"
  exit
fi
rm compiler.txt
#
gfortran file_transpose.o
if [ $? -ne 0 ]; then
  echo "Errors occurred while linking file_transpose.o"
  exit
fi
rm file_transpose.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/file_transpose
#
echo "Program installed as ~/bin/$ARCH/file_transpose"
