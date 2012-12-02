#!/bin/bash
#
gfortran -c args.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling args.f90"
  exit
fi
rm compiler.txt
#
F90 args.o
if [ $? -ne 0 ]; then
  echo "Errors loading args.o"
  exit
fi
rm args.o
#
mv a.out ~/bin/$ARCH/args
#
echo "Executable installed as ~/bin/$ARCH/args"
