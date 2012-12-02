#!/bin/bash
#
gfortran -c -g nco_tetrahedron_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nco_tetrahedron_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran nco_tetrahedron_prb.o -L$HOME/lib/$ARCH -lnco_tetrahedron
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nco_tetrahedron_prb.o"
  exit
fi
rm nco_tetrahedron_prb.o
#
mv a.out nco_tetrahedron_prb
./nco_tetrahedron_prb > nco_tetrahedron_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nco_tetrahedron_prb"
  exit
fi
rm nco_tetrahedron_prb
#
echo "Test program output written to nco_tetrahedron_prb_output.txt."
