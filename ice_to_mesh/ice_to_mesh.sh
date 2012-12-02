#!/bin/bash
#
g95 -c -g -I /usr/local/include ice_to_mesh.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling ice_to_mesh.f90"
  exit
fi
rm compiler.txt
#
g95 ice_to_mesh.o -L/usr/local/lib -lnetcdf
if [ $? -ne 0 ]; then
  echo "Errors while loading ice_to_mesh.o"
  exit
fi
rm ice_to_mesh.o
#
mv a.out ~/bin/$ARCH/ice_to_mesh
#
echo "Executable installed as ~/bin/$ARCH/ice_to_mesh"
