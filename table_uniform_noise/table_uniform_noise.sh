#!/bin/bash
#
gfortran -c -g table_uniform_noise.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_uniform_noise.f90"
  exit
fi
rm compiler.txt
#
gfortran table_uniform_noise.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_uniform_noise.o"
  exit
fi
rm table_uniform_noise.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_uniform_noise
#
echo "Executable installed as ~/bin/$ARCH/table_uniform_noise"
