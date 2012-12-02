#!/bin/bash
#
gfortran -c -g tri_surface_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tri_surface_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran tri_surface_io_prb.o -L$HOME/lib/$ARCH -ltri_surface_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tri_surface_io_prb.o"
  exit
fi
rm tri_surface_io_prb.o
#
mv a.out tri_surface_io_prb
./tri_surface_io_prb > tri_surface_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tri_surface_io_prb"
  exit
fi
rm tri_surface_io_prb
#
echo "Test program output written to tri_surface_io_prb_output.txt."
