#!/bin/bash
#
gfortran -c -g sphere_delaunay.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_delaunay.f90"
  exit
fi
rm compiler.txt
#
gfortran sphere_delaunay.o -L$HOME/lib/$ARCH -lstripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_delaunay.o"
  exit
fi
rm sphere_delaunay.o
#
mv a.out ~/bin/$ARCH/sphere_delaunay
#
echo "Executable installed as ~/bin/$ARCH/sphere_delaunay"
