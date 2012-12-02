#!/bin/bash
#
gfortran -c -g fd2d_predator_prey_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd2d_predator_prey_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran $HOME/lib/$ARCH/fd2d_predator_prey.o fd2d_predator_prey_prb.o -L$HOME/lib/$ARCH -ldlap
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd2d_predator_prey_prb.o"
  exit
fi
rm fd2d_predator_prey_prb.o
#
mv a.out fd2d_predator_prey_prb_gmres
./fd2d_predator_prey_prb_gmres < fd2d_predator_prey_prb_gmres_input.txt > fd2d_predator_prey_prb_gmres_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fd2d_predator_prey_prb_gmres"
  exit
fi
rm fd2d_predator_prey_prb_gmres
#
echo "Test program output written to fd2d_predator_prey_prb_gmres_output.txt."
