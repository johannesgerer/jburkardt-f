#!/bin/bash
#
gfortran -c -g table_orthonormalize.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_orthonormalize.f90"
  exit
fi
rm compiler.txt
gfortran table_orthonormalize.o -llapack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_orthonormalize.o"
  exit
fi
rm table_orthonormalize.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_orthonormalize
#
echo "Executable installed as ~/bin/$ARCH/table_orthonormalize"
