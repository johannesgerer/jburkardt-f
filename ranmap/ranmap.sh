#!/bin/bash
#
gfortran -c -g ranmap.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ranmap.f90"
  exit
fi
rm compiler.txt
#
gfortran ranmap.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ranmap.o"
  exit
fi
#
rm ranmap.o
mv a.out ~/bin/$ARCH/ranmap
#
echo "Executable installed as ~/bin/$ARCH/ranmap"
