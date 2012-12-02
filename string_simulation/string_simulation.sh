#!/bin/bash
#
gfortran -c string_simulation.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling string_simulation.f90"
  exit
fi
rm compiler.txt
#
gfortran string_simulation.o
if [ $? -ne 0 ]; then
  echo "Errors linking string_simulation.o"
  exit
fi
#
rm string_simulation.o
mv a.out ~/bin/$ARCH/string_simulation
#
echo "Executable installed as ~/bin/$ARCH/string_simulation"
