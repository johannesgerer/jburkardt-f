#!/bin/bash
#
gfortran -c -g sphere_grid_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_grid_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sphere_grid_prb.o -L$HOME/lib/$ARCH -lsphere_grid
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_grid_prb.o"
  exit
fi
rm sphere_grid_prb.o
#
mv a.out sphere_grid_prb
./sphere_grid_prb > sphere_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_grid_prb"
  exit
fi
rm sphere_grid_prb
#
echo "Test program output written to sphere_grid_prb_output.txt."
