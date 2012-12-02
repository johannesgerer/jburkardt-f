#!/bin/bash
#
gfortran -c -g ball_volume_monte_carlo.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_volume_monte_carlo.f90"
  exit
fi
rm compiler.txt
#
gfortran ball_volume_monte_carlo.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ball_volume_monte_carlo.o"
  exit
fi
rm ball_volume_monte_carlo.o
#
mv a.out ~/bin/$ARCH/ball_volume_monte_carlo
#
echo "Executable installed as ~/bin/$ARCH/ball_volume_monte_carlo"
