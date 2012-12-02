#!/bin/bash
#
gfortran -c -g int_exactness_hermite.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness_hermite.f90"
  exit
fi
rm compiler.txt
#
gfortran int_exactness_hermite.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness_hermite.o"
  exit
fi
rm int_exactness_hermite.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/int_exactness_hermite
#
echo "The program has been installed as ~/bin/$ARCH/int_exactness_hermite."
