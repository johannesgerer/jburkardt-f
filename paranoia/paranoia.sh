#!/bin/bash
#
gfortran -c paranoia.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling paranoia.f90"
  exit
fi
rm compiler.txt
#
gfortran paranoia.o
if [ $? -ne 0 ]; then
  echo "Errors loading paranoia.o"
  exit
fi
rm paranoia.o
#
mv a.out ~/bin/$ARCH/paranoia
#
echo "Executable installed as ~/bin/$ARCH/paranoia"
