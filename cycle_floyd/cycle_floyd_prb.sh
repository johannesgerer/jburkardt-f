#!/bin/bash
#
gfortran -c -g cycle_floyd_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cycle_floyd_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cycle_floyd_prb.o -L$HOME/lib/$ARCH -lcycle_floyd
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cycle_floyd_prb.o"
  exit
fi
rm cycle_floyd_prb.o
#
mv a.out cycle_floyd_prb
./cycle_floyd_prb > cycle_floyd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cycle_floyd_prb"
  exit
fi
rm cycle_floyd_prb
#
echo "Test program output written to cycle_floyd_prb_output.txt."
