#!/bin/bash
#
gfortran -c -g testpack.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling testpack.f90"
  exit
fi
rm compiler.txt
#
gfortran testpack.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading testpack.o"
  exit
fi
rm testpack.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/testpack
#
echo "Executable installed as ~/bin/$ARCH/testpack"
