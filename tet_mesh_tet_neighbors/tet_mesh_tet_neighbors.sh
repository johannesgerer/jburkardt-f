#!/bin/bash
#
F90 -c -g tet_mesh_tet_neighbors.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling tet_mesh_tet_neighbors.f90"
  exit
fi
rm compiler.txt
#
F90 tet_mesh_tet_neighbors.o
if [ $? -ne 0 ]; then
  echo "Error loading tet_mesh_tet_neighbors.o
  exit
fi
rm tet_mesh_tet_neighbors.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/tet_mesh_tet_neighbors
#
echo "Executable installed as ~/bin/$ARCH/tet_mesh_tet_neighbors"
