#!/bin/bash
#
gfortran -c -g rk4_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rk4_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran rk4_prb.o -L$HOME/lib/$ARCH -lrk4
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rk4_prb.o"
  exit
fi
rm rk4_prb.o
#
mv a.out rk4_prb
./rk4_prb > rk4_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rk4_prb"
  exit
fi
rm rk4_prb
#
echo "Test program output written to rk4_prb_output.txt."
