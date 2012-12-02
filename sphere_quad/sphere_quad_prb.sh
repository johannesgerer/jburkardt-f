#!/bin/bash
#
gfortran -c -g sphere_quad_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_quad_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sphere_quad_prb.o -L$HOME/lib/$ARCH -lsphere_quad
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_quad_prb.o"
  exit
fi
rm sphere_quad_prb.o
#
mv a.out sphere_quad_prb
./sphere_quad_prb > sphere_quad_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_quad_prb"
  exit
fi
rm sphere_quad_prb
#
echo "Test program output written to sphere_quad_prb_output.txt."
