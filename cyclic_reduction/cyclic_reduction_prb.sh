#!/bin/bash
#
gfortran -c -g cyclic_reduction_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cyclic_reduction_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cyclic_reduction_prb.o -L$HOME/lib/$ARCH -lcyclic_reduction
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cyclic_reduction_prb.o"
  exit
fi
rm cyclic_reduction_prb.o
#
mv a.out cyclic_reduction_prb
./cyclic_reduction_prb > cyclic_reduction_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cyclic_reduction_prb"
  exit
fi
rm cyclic_reduction_prb
#
echo "Test program output written to cyclic_reduction_prb_output.txt."
