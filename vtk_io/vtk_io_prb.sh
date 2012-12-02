#!/bin/bash
#
gfortran -c -g vtk_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vtk_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran vtk_io_prb.o -L$HOME/lib/$ARCH -lvtk_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vtk_io_prb.o"
  exit
fi
rm vtk_io_prb.o
#
mv a.out vtk_io_prb
./vtk_io_prb > vtk_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running vtk_io_prb"
  exit
fi
rm vtk_io_prb
#
echo "Test program output written to vtk_io_prb_output.txt."
