#!/bin/bash
#
gfortran -c md1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md1.f90"
  exit
fi
rm compiler.txt
#
gfortran md1.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md1.o"
  exit
fi
rm md1.o
#
mv a.out ~/bin/$ARCH/md1
#
echo "Executable installed as ~/bin/$ARCH/md1"
