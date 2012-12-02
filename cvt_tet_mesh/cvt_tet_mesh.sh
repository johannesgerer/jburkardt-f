#!/bin/bash
#
gfortran -c -g cvt_tet_mesh.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_tet_mesh.csh"
  exit
fi
rm compiler.txt
#
gfortran cvt_tet_mesh.o -L$HOME/lib/$ARCH -ltest_tet_mesh
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_tet_mesh.o"
  exit
fi
rm cvt_tet_mesh.o
#
mv a.out cvt_tet_mesh
./cvt_tet_mesh > cvt_tet_mesh_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_tet_mesh"
  exit
fi
rm cvt_tet_mesh
#
echo "Program output written to cvt_tet_mesh_output.txt"
