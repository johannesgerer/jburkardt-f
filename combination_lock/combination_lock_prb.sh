#!/bin/bash
#
gfortran -c -g combination_lock_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling combination_lock_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran combination_lock_prb.o -L$HOME/lib/$ARCH -lcombination_lock
if [ $? -ne 0 ]; then
  echo "Errors linking and loading combination_lock_prb.o"
  exit
fi
rm combination_lock_prb.o
#
mv a.out combination_lock_prb
./combination_lock_prb > combination_lock_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running combination_lock_prb"
  exit
fi
rm combination_lock_prb
#
echo "Program output written to combination_lock_prb_output.txt"
