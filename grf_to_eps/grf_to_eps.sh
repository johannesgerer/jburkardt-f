#!/bin/bash
#
gfortran -c -g grf_to_eps.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grf_to_eps.f90"
  exit
fi
rm compiler.txt
#
gfortran grf_to_eps.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grf_to_eps.o"
  exit
fi
#
rm grf_to_eps.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/grf_to_eps
#
echo "Program installed as ~/bin/$ARCH/grf_to_eps"
