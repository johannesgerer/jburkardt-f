#!/bin/bash
#
gfortran -c -g ball_volume_quad.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_volume_quad.f90"
  exit
fi
rm compiler.txt
#
gfortran ball_volume_quad.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ball_volume_quad.o"
  exit
fi
rm ball_volume_quad.o
#
mv a.out ~/bin/$ARCH/ball_volume_quad
#
echo "Executable installed as ~/bin/$ARCH/ball_volume_quad"
