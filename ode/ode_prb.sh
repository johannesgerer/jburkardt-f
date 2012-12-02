#!/bin/bash
#
gfortran -c -g ode_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ode_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ode_prb.o -L$HOME/lib/$ARCH -lode
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ode_prb.o"
  exit
fi
rm ode_prb.o
#
mv a.out ode_prb
./ode_prb > ode_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ode_prb"
  exit
fi
rm ode_prb
#
echo "Test program output written to ode_prb_output.txt."
