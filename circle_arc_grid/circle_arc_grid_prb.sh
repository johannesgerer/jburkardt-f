#!/bin/bash
#
gfortran -c -g circle_arc_grid_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_arc_grid_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran circle_arc_grid_prb.o -L$HOME/lib/$ARCH -lcircle_arc_grid
if [ $? -ne 0 ]; then
  echo "Errors linking and loading circle_arc_grid_prb.o"
  exit
fi
rm circle_arc_grid_prb.o
#
mv a.out circle_arc_grid_prb
./circle_arc_grid_prb > circle_arc_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running circle_arc_grid_prb"
  exit
fi
rm circle_arc_grid_prb
#
echo "Program output written to circle_arc_grid_prb_output.txt"
