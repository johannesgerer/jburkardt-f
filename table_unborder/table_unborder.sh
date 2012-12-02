#!/bin/bash
#
gfortran -c -g table_unborder.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_unborder.f90"
  exit
fi
rm compiler.txt
#
gfortran table_unborder.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_unborder.o"
  exit
fi
rm table_unborder.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_unborder
#
echo "Executable installed as ~/bin/$ARCH/table_unborder"
