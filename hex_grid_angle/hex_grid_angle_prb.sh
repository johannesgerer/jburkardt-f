#!/bin/bash
#
gfortran -c -g hex_grid_angle_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hex_grid_angle_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran hex_grid_angle_prb.o -L$HOME/lib/$ARCH -lhex_grid_angle
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hex_grid_angle_prb.o"
  exit
fi
rm hex_grid_angle_prb.o
#
mv a.out hex_grid_angle_prb
./hex_grid_angle_prb > hex_grid_angle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hex_grid_angle_prb"
  exit
fi
rm hex_grid_angle_prb
#
echo "Test program output written to hex_grid_angle_prb_output.txt."
