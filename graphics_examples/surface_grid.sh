#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran surface_grid.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling surface_grid.f90"
  exit
fi
rm compiler.txt
#
gfortran surface_grid.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading surface_grid.o"
  exit
fi
rm surface_grid.o
#
mv a.out surface_grid
./surface_grid > surface_grid_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running surface_grid"
  exit
fi
rm surface_grid
#
echo "Program output written to surface_grid_output.txt"
