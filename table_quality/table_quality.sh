#!/bin/bash
#
gfortran -c -g table_quality.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling table_quality.f90"
  exit
fi
rm compiler.txt
#
gfortran table_quality.o
if [ $? -ne 0 ]; then
  echo "Error loading table_quality.o"
  exit
fi
rm table_quality.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_quality
#
echo "Executable installed as ~/bin/$ARCH/table_quality"
