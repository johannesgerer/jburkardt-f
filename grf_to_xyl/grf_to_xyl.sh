#!/bin/bash
#
gfortran -c -g grf_to_xyl.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grf_to_xyl.f90"
  exit
fi
rm compiler.txt
#
gfortran grf_to_xyl.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grf_to_xyl.o"
  exit
fi
#
rm grf_to_xyl.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/grf_to_xyl
#
echo "Program installed as ~/bin/$ARCH/grf_to_xyl"
