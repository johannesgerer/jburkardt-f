#!/bin/bash
#
gfortran -c -g hex_grid_triangulate.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hex_grid_triangulate.f90"
  exit
fi
rm compiler.txt
#
gfortran hex_grid_triangulate.o -L$HOME/lib/$ARCH -lhex_grid -ltest_triangulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hex_grid_triangulate.o"
  exit
fi
rm hex_grid_triangulate.o
#
mv a.out hex_grid_triangulate
./hex_grid_triangulate > hex_grid_triangulate_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hex_grid_triangulate"
  exit
fi
rm hex_grid_triangulate
#
echo "Program output written to hex_grid_triangulate_output.txt"
