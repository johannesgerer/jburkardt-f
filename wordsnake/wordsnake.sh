#!/bin/bash
#
gfortran -c wordsnake.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wordsnake.f90"
  exit
fi
rm compiler.txt
#
gfortran wordsnake.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wordsnake.o"
  exit
fi
rm wordsnake.o
#
mv a.out ~/bin/$ARCH/wordsnake
#
echo "Executable installed as ~/bin/$ARCH/wordsnake"
