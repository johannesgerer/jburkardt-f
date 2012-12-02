#!/bin/bash
#
gfortran -c -g table_stats.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_stats.f90"
  exit
fi
rm compiler.txt
#
gfortran table_stats.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_stats.o"
  exit
fi
rm table_stats.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_stats
#
echo "Executable installed as ~/bin/$ARCH/table_stats"
