#!/bin/bash
#
gfortran -c -g table_border.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_border.f90"
  exit
fi
rm compiler.txt
#
gfortran table_border.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_border.o"
  exit
fi
rm table_border.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_border
#
echo "Exectuable installed as ~/bin/$ARCH/table_border"
