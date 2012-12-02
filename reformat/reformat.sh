#!/bin/bash
#
gfortran -c -g reformat.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling reformat.f90"
  exit
fi
rm compiler.txt
#
gfortran reformat.o
if [ $? -ne 0 ]; then
  echo "Errors loading reformat.o"
  exit
fi
rm reformat.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/reformat
#
echo "Executable installed as ~/bin/$ARCH/reformat"
