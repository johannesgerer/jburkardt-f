#!/bin/bash
#
gfortran -c -g ising_3d_simulation.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ising_3d_simulation.f90"
  exit
fi
rm compiler.txt
#
gfortran ising_3d_simulation.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ising_3d_simulation.o"
  exit
fi
rm ising_3d_simulation.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/ising_3d_simulation
#
echo "The program has been installed as ~/bin/$ARCH/ising_3d_simulation."
