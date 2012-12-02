#!/bin/bash
#
gfortran -c md.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md.f90"
  exit
fi
rm compiler.txt
#
gfortran md.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md.o"
  exit
fi
rm md.o
#
mv a.out ~/bin/$ARCH/md
#
echo "Executable installed as ~/bin/$ARCH/md"
