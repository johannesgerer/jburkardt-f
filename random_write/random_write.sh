#!/bin/bash
#
gfortran -c random_write.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_write.f90"
  exit
fi
rm compiler.txt
#
gfortran random_write.o
if [ $? -ne 0 ]; then
  echo "Errors loading random_write.o"
  exit
fi
rm random_write.o
#
mv a.out ~/bin/$ARCH/random_write
#
echo "Executable installed as ~/bin/$ARCH/random_write"
