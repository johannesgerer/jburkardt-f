#!/bin/bash
#
gfortran -c -g int_exactness_gen_hermite.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling int_exactness_gen_hermite.f90"
  exit
fi
rm compiler.txt
#
gfortran int_exactness_gen_hermite.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading int_exactness_gen_hermite.o"
  exit
fi
rm int_exactness_gen_hermite.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/int_exactness_gen_hermite
#
echo "Executable installed as ~/bin/$ARCH/int_exactness_gen_hermite"
