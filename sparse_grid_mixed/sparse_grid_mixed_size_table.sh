#!/bin/bash
#
gfortran -c -g sparse_grid_mixed_size_table.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_mixed_size_table.f90"
  exit
fi
rm compiler.txt
#
gfortran sparse_grid_mixed_size_table.o -L$HOME/lib/$ARCH -lsparse_grid_mixed -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_mixed_size_table.o"
  exit
fi
rm sparse_grid_mixed_size_table.o
#
mv a.out sparse_grid_mixed_size_table
./sparse_grid_mixed_size_table > sparse_grid_mixed_size_table_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_grid_mixed_size_table"
  exit
fi
rm sparse_grid_mixed_size_table
#
echo "Program output written to sparse_grid_mixed_size_table_output.txt"
