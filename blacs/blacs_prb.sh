#!/bin/bash
#
mpif90 -c blacs_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling blacs_prb.f90"
  exit
fi
rm compiler.txt
#
mpif90 blacs_prb.o -L ~/lib/$ARCH -lblacs -lmpi
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading blacs_prb.o"
  exit
fi
rm blacs_prb.o
#
mv a.out blacs_prb
mpirun -v -np 4 ./blacs_prb > blacs_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running blacs_prb"
  exit
fi
rm blacs_prb
#
echo "Program output written to blacs_prb_output.txt"
