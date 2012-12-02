#!/bin/bash
#
gfortran -c -g extract.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling extract.f90"
  exit
fi
rm compiler.txt
#
gfortran extract.o
if [ $? -ne 0 ]; then
  echo "Errors loading extract.o"
  exit
fi
rm extract.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/extract
#
echo "The program was installed as ~/bin/$ARCH/extract"
