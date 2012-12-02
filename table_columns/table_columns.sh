#!/bin/bash
#
gfortran -c -g table_columns.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_columns.f90"
  exit
fi
rm compiler.txt
#
gfortran table_columns.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_columns.o"
  exit
fi
rm table_columns.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_columns
#
echo "Program installed as ~/bin/$ARCH/table_columns"
