#!/bin/bash
#
gfortran -c -g stripper.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stripper.f90"
  exit
fi
rm compiler.txt
#
gfortran stripper.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stripper.o"
  exit
fi
rm stripper.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/stripper
#
echo "Executable installed as ~/bin/$ARCH/stripper"
