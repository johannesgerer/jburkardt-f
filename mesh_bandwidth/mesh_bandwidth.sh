#!/bin/bash
#
gfortran -c -g mesh_bandwidth.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mesh_bandwidth.f90"
  exit
fi
rm compiler.txt
#
gfortran mesh_bandwidth.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mesh_bandwidth.o"
  exit
fi
rm mesh_bandwidth.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/mesh_bandwidth
#
echo "Executable installed as ~/bin/$ARCH/mesh_bandwidth"
