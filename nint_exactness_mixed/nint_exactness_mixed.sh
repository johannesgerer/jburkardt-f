#!/bin/bash
#
gfortran -c -g nint_exactness_mixed.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nint_exactness_mixed.f90"
  exit
fi
rm compiler.txt
#
gfortran nint_exactness_mixed.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nint_exactness_mixed.o"
  exit
fi
rm nint_exactness_mixed.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/nint_exactness_mixed
#
echo "Executable installed as ~/bin/$ARCH/nint_exactness_mixed"
