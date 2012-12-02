#!/bin/bash
#
gfortran -c -g int_exactness_gegenbauer.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness_gegenbauer.f90"
  exit
fi
rm compiler.txt
#
gfortran int_exactness_gegenbauer.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness_gegenbauer.o"
  exit
fi
rm int_exactness_gegenbauer.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/int_exactness_gegenbauer
#
echo "Executable installed as ~/bin/$ARCH/int_exactness_gegenbauer"
