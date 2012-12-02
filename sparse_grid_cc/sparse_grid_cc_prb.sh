#!/bin/bash
#
gfortran -c -g sparse_grid_cc_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_cc_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sparse_grid_cc_prb.o -L$HOME/lib/$ARCH -lsparse_grid_cc
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_cc_prb.o"
  exit
fi
rm sparse_grid_cc_prb.o
#
mv a.out sparse_grid_cc_prb
./sparse_grid_cc_prb > sparse_grid_cc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_grid_cc_prb"
  exit
fi
rm sparse_grid_cc_prb
#
echo "Test program output written to sparse_grid_cc_prb_output.txt."
