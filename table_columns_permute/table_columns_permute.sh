#!/bin/bash
#
gfortran -c -g table_columns_permute.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_columns_permute.f90"
  exit
fi
rm compiler.txt
#
gfortran table_columns_permute.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_columns_permute.o"
  exit
fi
rm table_columns_permute.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_columns_permute
#
echo "Executable installed as ~/bin/$ARCH/table_columns_permute"
