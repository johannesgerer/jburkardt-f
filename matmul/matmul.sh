#!/bin/bash
#
gfortran -c -g matmul.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matmul.f90"
  exit
fi
rm compiler.txt
#
gfortran matmul.o -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matmul.o"
  exit
fi
rm matmul.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/matmul
#
echo "Executable installed as ~/bin/$ARCH/matmul"
