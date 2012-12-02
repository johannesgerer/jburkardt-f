#!/bin/bash
#
gfortran -c -g fd1d_predator_prey_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_predator_prey_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran fd1d_predator_prey_prb.o -L$HOME/lib/$ARCH -lfd1d_predator_prey
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd1d_predator_prey_prb.o"
  exit
fi
rm fd1d_predator_prey_prb.o
#
mv a.out fd1d_predator_prey_prb
./fd1d_predator_prey_prb < fd1d_predator_prey_input.txt > fd1d_predator_prey_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fd1d_predator_prey_prb"
  exit
fi
rm fd1d_predator_prey_prb
#
echo "Test program output written to fd1d_predator_prey_prb_output.txt."
