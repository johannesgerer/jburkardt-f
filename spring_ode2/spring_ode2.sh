#!/bin/bash
#
gfortran -c spring_ode2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spring_ode2.f90"
  exit
fi
rm compiler.txt
#
gfortran spring_ode2.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spring_ode2.o"
  exit
fi
rm spring_ode2.o
#
mv a.out ~/bin/$ARCH/spring_ode2
#
echo "Executable installed as ~/bin/$ARCH/spring_ode2"
