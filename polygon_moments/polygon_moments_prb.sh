#!/bin/bash
#
gfortran -c -g polygon_moments_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_moments_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran polygon_moments_prb.o -L$HOME/lib/$ARCH -lpolygon_moments
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygon_moments_prb.o"
  exit
fi
rm polygon_moments_prb.o
#
mv a.out polygon_moments_prb
./polygon_moments_prb > polygon_moments_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polygon_moments_prb"
  exit
fi
rm polygon_moments_prb
#
echo "Test program output written to polygon_moments_prb_output.txt."
