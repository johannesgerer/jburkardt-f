#!/bin/bash
#
gfortran -c -g r16_int_exactness_gen_hermite.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r16_int_exactness_gen_hermite.f90"
  exit
fi
rm compiler.txt
#
gfortran r16_int_exactness_gen_hermite.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading r16_int_exactness_gen_hermite.o"
  exit
fi
rm r16_int_exactness_gen_hermite.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/r16_int_exactness_gen_hermite
#
echo "Executable installed as ~/bin/$ARCH/r16_int_exactness_gen_hermite"
