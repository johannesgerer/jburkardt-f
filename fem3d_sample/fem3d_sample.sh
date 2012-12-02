#!/bin/bash
#
gfortran -c fem3d_sample.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem3d_sample.f90"
  exit
fi
rm compiler.txt
#
gfortran fem3d_sample.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem3d_sample.o"
  exit
fi
#
rm fem3d_sample.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/fem3d_sample
#
echo "Program installed as ~/bin/$ARCH/fem3d_sample"
