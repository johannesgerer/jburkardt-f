#!/bin/bash
#
gfortran -c spring_ode.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spring_ode.f90"
  exit
fi
rm compiler.txt
#
gfortran spring_ode.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spring_ode.o"
  exit
fi
rm spring_ode.o
#
mv a.out ~/bin/$ARCH/spring_ode
#
echo "Executable installed as ~/bin/$ARCH/spring_ode"
