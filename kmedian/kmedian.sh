#!/bin/bash
#
gfortran -c -g kmedian.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kmedian.f90"
  exit
fi
rm compiler.txt
#
gfortran kmedian.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading kmedian.o"
  exit
fi
rm kmedian.o
#
mv a.out ~/bin/$ARCH/kmedian
#
echo "Executable installed as ~/bin/$ARCH/kmedian"
