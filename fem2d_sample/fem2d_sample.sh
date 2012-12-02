#!/bin/bash
#
gfortran -c fem2d_sample.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_sample.f90"
  exit
fi
rm compiler.txt
#
gfortran fem2d_sample.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_sample.o"
  exit
fi
#
rm fem2d_sample.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/fem2d_sample
#
echo "Program installed as ~/bin/$ARCH/fem2d_sample"
