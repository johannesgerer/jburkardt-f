#!/bin/bash
#
gfortran -c -g triangulation_triangle_neighbors.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_triangle_neighbors.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_triangle_neighbors.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_triangle_neighbors.o"
  exit
fi
rm triangulation_triangle_neighbors.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_triangle_neighbors
#
echo "Executable installed as ~/bin/$ARCH/triangulation_triangle_neighbors"
