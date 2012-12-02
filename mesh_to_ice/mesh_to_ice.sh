#!/bin/bash
#
g95 -c -g -I /usr/local/include mesh_to_ice.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling mesh_to_ice.f90"
  exit
fi
rm compiler.txt
#
g95 mesh_to_ice.o -L/usr/local/lib -lnetcdf
if [ $? -ne 0 ]; then
  echo "Errors while loading mesh_to_ice.o"
  exit
fi
rm mesh_to_ice.o
#
mv a.out ~/bin/$ARCH/mesh_to_ice
#
echo "Executable installed as ~/bin/$ARCH/mesh_to_ice"
