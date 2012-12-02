#!/bin/bash
#
gfortran -c -g linpack_z_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_z_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran linpack_z_prb.o -L$HOME/lib/$ARCH -llinpack_z -lblas1_z -lblas1_d
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_z_prb.o"
  exit
fi
rm linpack_z_prb.o
#
mv a.out linpack_z_prb
./linpack_z_prb > linpack_z_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linpack_z_prb"
  exit
fi
rm linpack_z_prb
#
echo "Test program output written to linpack_z_prb_output.txt."
