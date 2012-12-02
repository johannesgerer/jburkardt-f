#!/bin/bash
#
gfortran -c -g quad_mesh_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_mesh_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran quad_mesh_prb.o -L$HOME/lib/$ARCH -lquad_mesh
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quad_mesh_prb.o"
  exit
fi
rm quad_mesh_prb.o
#
mv a.out quad_mesh_prb
./quad_mesh_prb > quad_mesh_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quad_mesh_prb"
  exit
fi
rm quad_mesh_prb
#
echo "Test program output written to quad_mesh_prb_output.txt."
