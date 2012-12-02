#!/bin/bash
#
gfortran -c -g table_merge.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling merge.f90"
  exit
fi
rm compiler.txt
#
gfortran table_merge.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_merge.o"
  exit
fi
rm table_merge.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_merge
#
echo "Program installed as ~/bin/$ARCH/table_merge"
