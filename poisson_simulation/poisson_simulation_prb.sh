#!/bin/bash
#
gfortran -c -g poisson_simulation_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_simulation_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran poisson_simulation_prb.o -L$HOME/lib/$ARCH -lpoisson_simulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading poisson_simulation_prb.o"
  exit
fi
rm poisson_simulation_prb.o
#
mv a.out poisson_simulation_prb
./poisson_simulation_prb > poisson_simulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poisson_simulation_prb"
  exit
fi
rm poisson_simulation_prb
#
echo "Test program output written to poisson_simulation_prb_output.txt."
