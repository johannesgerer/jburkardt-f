#!/bin/bash
#
gfortran -c -g pce_ode_hermite_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pce_ode_hermite_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran pce_ode_hermite_prb.o -L$HOME/lib/$ARCH -lpce_ode_hermite
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pce_ode_hermite_prb.o"
  exit
fi
rm pce_ode_hermite_prb.o
#
mv a.out pce_ode_hermite_prb
./pce_ode_hermite_prb > pce_ode_hermite_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pce_ode_hermite_prb"
  exit
fi
rm pce_ode_hermite_prb
#
echo "Test program output written to pce_ode_hermite_prb_output.txt."
