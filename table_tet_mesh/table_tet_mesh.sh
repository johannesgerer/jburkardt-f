#!/bin/bash
#
gfortran -c -g table_tet_mesh.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_tet_mesh.f90"
  exit
fi
rm compiler.txt
#
gfortran table_tet_mesh.o -L$HOME/lib/$ARCH
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_tet_mesh.o"
  exit
fi
#
rm table_tet_mesh.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_tet_mesh
#
echo "Executable installed as ~/bin/$ARCH/table_tet_mesh"
