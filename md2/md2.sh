#!/bin/bash
#
gfortran -c md2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md2.f90"
  exit
fi
rm compiler.txt
#
gfortran md2.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md2.o"
  exit
fi
rm md2.o
#
mv a.out ~/bin/$ARCH/md2
#
echo "Executable installed as ~/bin/$ARCH/md2"
