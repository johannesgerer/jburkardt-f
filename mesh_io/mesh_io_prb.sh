#!/bin/bash
#
gfortran -c -g mesh_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mesh_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran mesh_io_prb.o -L$HOME/lib/$ARCH -lmesh_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mesh_io_prb.o"
  exit
fi
rm mesh_io_prb.o
#
mv a.out mesh_io_prb
./mesh_io_prb > mesh_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mesh_io_prb"
  exit
fi
rm mesh_io_prb
#
echo "Test program output written to mesh_io_prb_output.txt."
