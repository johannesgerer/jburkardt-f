#!/bin/bash
#
gfortran -c -g inout.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling inout.f90"
  exit
fi
rm compiler.txt
#
gfortran ~/lib/$ARCH/fem2d_stokes.o inout.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_stokes.o + inout.o"
  exit
fi
rm inout.o
#
chmod ugo+x a.out
mv a.out inout
./inout nodes6.txt triangles6.txt > inout_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running inout."
  exit
fi
rm inout
#
echo "Program output written to inout_output.txt"
