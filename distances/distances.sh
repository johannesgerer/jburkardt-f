#!/bin/bash
#
gfortran -c distances.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling distances.f90"
  exit
fi
rm compiler.txt
#
gfortran distances.o
if [ $? -ne 0 ]; then
  echo "Errors loading distances.o"
  exit
fi
#
rm distances.o
mv a.out ~/bin/$ARCH/distances
#
echo "Executable installed as ~/bin/$ARCH/distances"
