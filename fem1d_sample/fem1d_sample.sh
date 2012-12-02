#!/bin/bash
#
gfortran -c fem1d_sample.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_sample.f90"
  exit
fi
rm compiler.txt
#
gfortran fem1d_sample.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_sample.o"
  exit
fi
#
rm fem1d_sample.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/fem1d_sample
#
echo "Program installed as ~/bin/$ARCH/fem1d_sample"
