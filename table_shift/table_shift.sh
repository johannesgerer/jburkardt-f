#!/bin/bash
#
gfortran -c -g table_shift.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_shift.f90"
  exit
fi
rm compiler.txt
#
gfortran table_shift.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_shift.o"
  exit
fi
rm table_shift.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_shift
#
echo "The executable has been installed as ~/bin/$ARCH/table_shift"
