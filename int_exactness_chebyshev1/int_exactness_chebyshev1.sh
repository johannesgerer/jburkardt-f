#!/bin/bash
#
gfortran -c -g int_exactness_chebyshev1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness_chebyshev1.f90"
  exit
fi
rm compiler.txt
#
gfortran int_exactness_chebyshev1.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness_chebyshev1.o"
  exit
fi
rm int_exactness_chebyshev1.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/int_exactness_chebyshev1
#
echo "Executable installed as ~/bin/$ARCH/int_exactness_chebyshev1"
