#!/bin/bash
#
gfortran -c -g table_scale.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_scale.f90"
  exit
fi
rm compiler.txt
#
gfortran table_scale.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_scale.o"
  exit
fi
rm table_scale.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_scale
#
echo "Executable installed as ~/bin/$ARCH/table_scale"
