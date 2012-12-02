#!/bin/bash
#
gfortran -c -g geometry_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geometry_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran geometry_prb.o -L$HOME/lib/$ARCH -lgeometry
if [ $? -ne 0 ]; then
  echo "Errors linking and loading geometry_prb.o"
  exit
fi
rm geometry_prb.o
#
mv a.out geometry_prb
./geometry_prb > geometry_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running geometry_prb"
  exit
fi
rm geometry_prb
#
echo "Test program output written to geometry_prb_output.txt."
