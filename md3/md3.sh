#!/bin/bash
#
gfortran -c md3.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md3.f90"
  exit
fi
rm compiler.txt
#
gfortran md3.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md3.o"
  exit
fi
rm md3.o
#
mv a.out ~/bin/$ARCH/md3
#
echo "Executable installed as ~/bin/$ARCH/md3"
