#!/bin/bash
#
gfortran -c -g quad_mesh_rcm.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_mesh_rcm.f90"
  exit
fi
rm compiler.txt
#
gfortran quad_mesh_rcm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quad_mesh_rcm.o"
  exit
fi
rm quad_mesh_rcm.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/quad_mesh_rcm
#
echo "Executable installed as ~/bin/$ARCH/quad_mesh_rcm"
