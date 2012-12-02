#!/bin/bash
#
gfortran -c -g pwl_interp_2d_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran pwl_interp_2d_prb.o -L$HOME/lib/$ARCH -lpwl_interp_2d -ltest_interp_2d -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pwl_interp_2d_prb.o"
  exit
fi
rm pwl_interp_2d_prb.o
#
mv a.out pwl_interp_2d_prb
./pwl_interp_2d_prb > pwl_interp_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pwl_interp_2d_prb"
  exit
fi
rm pwl_interp_2d_prb
#
echo "Test program output written to pwl_interp_2d_prb_output.txt."
