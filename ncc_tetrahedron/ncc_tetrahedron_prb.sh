#!/bin/bash
#
gfortran -c -g ncc_tetrahedron_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ncc_tetrahedron_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ncc_tetrahedron_prb.o -L$HOME/lib/$ARCH -lncc_tetrahedron
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ncc_tetrahedron_prb.o"
  exit
fi
rm ncc_tetrahedron_prb.o
#
mv a.out ncc_tetrahedron_prb
./ncc_tetrahedron_prb > ncc_tetrahedron_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ncc_tetrahedron_prb"
  exit
fi
rm ncc_tetrahedron_prb
#
echo "Test program output written to ncc_tetrahedron_prb_output.txt."
