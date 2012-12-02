#!/bin/bash
#
gfortran -c -g sphere_voronoi.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_voronoi.f90"
  exit
fi
rm compiler.txt
#
gfortran sphere_voronoi.o -L$HOME/lib/$ARCH -lstripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_voronoi.o"
  exit
fi
rm sphere_voronoi.o
#
mv a.out ~/bin/$ARCH/sphere_voronoi
#
echo "Executable installed as ~/bin/$ARCH/sphere_voronoi"
