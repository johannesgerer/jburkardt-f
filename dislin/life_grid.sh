#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran life_grid.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling life_grid.f90"
  exit
fi
rm compiler.txt
#
gfortran life_grid.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading life_grid.o"
  exit
fi
rm life_grid.o
#
mv a.out life_grid
./life_grid > life_grid_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running life_grid"
  exit
fi
rm life_grid
#
echo "Program output written to life_grid_output.txt"
