#!/bin/bash
#
gfortran -c -g table_voronoi.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_voronoi.f90"
  exit
fi
rm compiler.txt
#
gfortran table_voronoi.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_voronoi.o"
  exit
fi
rm table_voronoi.o
#
chmod ugo+x a.out
mv a.out $HOME/bin/$ARCH/table_voronoi
#
echo "Executable installed as $HOME/bin/$ARCH/table_voronoi"
