#!/bin/bash
#
gfortran -c -g sparse_grid_mixed_write_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_mixed_write_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sparse_grid_mixed_write_prb.o -L$HOME/lib/$ARCH -lsparse_grid_mixed -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_mixed_write_prb.o"
  exit
fi
rm sparse_grid_mixed_write_prb.o
#
mv a.out sparse_grid_mixed_write_prb
./sparse_grid_mixed_write_prb > sparse_grid_mixed_write_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_grid_mixed_write_prb"
  exit
fi
rm sparse_grid_mixed_write_prb
#
echo "Program output written to sparse_grid_mixed_write_prb_output.txt"
#
#  Move sparse grid files to dataset directory.
#
mv *_a.txt ../../datasets/sparse_grid_mixed
mv *_b.txt ../../datasets/sparse_grid_mixed
mv *_r.txt ../../datasets/sparse_grid_mixed
mv *_w.txt ../../datasets/sparse_grid_mixed
mv *_x.txt ../../datasets/sparse_grid_mixed
#
echo "Program output files moved to ../../datasets/sparse_grid_mixed"
