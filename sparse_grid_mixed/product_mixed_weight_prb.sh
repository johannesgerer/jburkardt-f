#!/bin/bash
#
gfortran -c -g product_mixed_weight_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling product_mixed_weight_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran product_mixed_weight_prb.o -L$HOME/lib/$ARCH -lsparse_grid_mixed -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading product_mixed_weight_prb.o"
  exit
fi
rm product_mixed_weight_prb.o
#
mv a.out product_mixed_weight_prb
./product_mixed_weight_prb > product_mixed_weight_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running product_mixed_weight_prb"
  exit
fi
rm product_mixed_weight_prb
#
echo "Program output written to product_mixed_weight_prb_output.txt"
