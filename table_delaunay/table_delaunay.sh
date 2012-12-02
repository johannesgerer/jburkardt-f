#!/bin/bash
#
gfortran -c -g table_delaunay.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_delaunay.f90"
  exit
fi
rm compiler.txt
#
gfortran table_delaunay.o -L$HOME/lib/$ARCH
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_delaunay.o"
  exit
fi
#
rm table_delaunay.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_delaunay
#
echo "Executable installed as ~/bin/$ARCH/table_delaunay"
