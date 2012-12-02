#!/bin/bash
#
gfortran -c -g table_barplot_ppma.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_barplot_ppma.f90"
  exit
fi
rm compiler.txt
#
gfortran table_barplot_ppma.o
if [ $? -ne 0 ]; then
  echo "Errors loading table_barplot_ppma.o"
  exit
fi
#
rm table_barplot_ppma.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_barplot_ppma
#
echo "Executable installed as ~/bin/$ARCH/table_barplot_ppma"
