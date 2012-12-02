#!/bin/bash
#
gfortran -c -g comp_next_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling comp_next_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran comp_next_prb.o -L$HOME/lib/$ARCH -lsparse_grid_mixed -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading comp_next_prb.o"
  exit
fi
rm comp_next_prb.o
#
mv a.out comp_next_prb
./comp_next_prb > comp_next_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running comp_next_prb"
  exit
fi
rm comp_next_prb
#
echo "Program output written to comp_next_prb_output.txt"
