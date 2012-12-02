#!/bin/bash
#
gfortran -c -g delaunay_lmap_2d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling delaunay_lmap_2d.f90"
  exit
fi
rm compiler.txt
#
gfortran delaunay_lmap_2d.o -L$HOME/lib/$ARCH
if [ $? -ne 0 ]; then
  echo "Errors linking and loading delaunay_lmap_2d.o"
  exit
fi
#
rm delaunay_lmap_2d.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/delaunay_lmap_2d
#
echo "Program installed as ~/bin/$ARCH/delaunay_lmap_2d"
