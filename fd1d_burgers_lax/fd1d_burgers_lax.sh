#!/bin/bash
#
gfortran -c -g fd1d_burgers_lax.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling fd1d_burgers_lax.f90"
  exit
fi
rm compiler.txt
#
gfortran fd1d_burgers_lax.o
if [ $? -ne 0 ]; then
  echo "Errors while loading fd1d_burgers_lax.o"
  exit
fi
rm fd1d_burgers_lax.o
#
mv a.out ~/bin/$ARCH/fd1d_burgers_lax
#
echo "Program installed as ~/bin/$ARCH/fd1d_burgers_lax"
