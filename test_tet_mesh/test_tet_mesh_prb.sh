#!/bin/bash
#
gfortran -c -g test_tet_mesh_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_tet_mesh_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_tet_mesh_prb.o -L$HOME/lib/$ARCH -ltest_tet_mesh
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_tet_mesh_prb.o"
  exit
fi
rm test_tet_mesh_prb.o
#
mv a.out test_tet_mesh_prb
./test_tet_mesh_prb > test_tet_mesh_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_tet_mesh_prb"
  exit
fi
rm test_tet_mesh_prb
#
echo "Test program output written to test_tet_mesh_prb_output.txt."
