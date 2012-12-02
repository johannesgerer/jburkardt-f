#!/bin/bash
#
gfortran -c -g sde_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sde_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sde_prb.o -L$HOME/lib/$ARCH -lsde -lqr_solve
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sde_prb.o"
  exit
fi
rm sde_prb.o
#
mv a.out sde_prb
./sde_prb > sde_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sde_prb"
  exit
fi
rm sde_prb
#
echo "Test program output written to sde_prb_output.txt."
