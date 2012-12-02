#!/bin/bash
#
gfortran -c -g int_exactness.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness.f90"
  exit
fi
rm compiler.txt
#
gfortran int_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness.o"
  exit
fi
rm int_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/int_exactness
#
echo "Executable installed as ~/bin/$ARCH/int_exactness"
