#!/bin/bash
#
gfortran -c -g ieee_uniform_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ieee_uniform_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ieee_uniform_prb.o -L$HOME/lib/$ARCH -lieee_uniform
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ieee_uniform_prb.o"
  exit
fi
rm ieee_uniform_prb.o
#
mv a.out ieee_uniform_prb
./ieee_uniform_prb > ieee_uniform_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ieee_uniform_prb"
  exit
fi
rm ieee_uniform_prb
#
echo "Test program output written to ieee_uniform_prb_output.txt."
