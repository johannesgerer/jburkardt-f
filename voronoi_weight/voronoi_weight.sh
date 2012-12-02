#!/bin/bash
#
gfortran -c -g voronoi_weight.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling voronoi_weight.f90"
  exit
fi
rm compiler.txt
#
gfortran voronoi_weight.o
if [ $? -ne 0 ]; then
  echo "Errors loading voronoi_weight.o"
  exit
fi
rm voronoi_weight.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/voronoi_weight
#
echo "Executable installed as ~/bin/$ARCH/voronoi_weight"
