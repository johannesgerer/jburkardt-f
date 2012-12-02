#!/bin/bash
#
gfortran -c -g fd2d_predator_prey.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd2d_predator_prey.f90"
  exit
fi
rm compiler.txt
#
mv fd2d_predator_prey.o ~/lib/$ARCH
#
echo "Object code installed as ~/lib/$ARCH/fd2d_predator_prey.o"
