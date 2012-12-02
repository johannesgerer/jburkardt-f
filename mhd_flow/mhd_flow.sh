#!/bin/bash
#
gfortran -c -g mhd_flow.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mhd_flow.f90"
  exit
fi
rm compiler.txt
#
gfortran mhd_flow.o -L$HOME/lib/$ARCH -llinpack_d -lblas1_d >& loader.txt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mhd_flow.o"
  exit
fi
rm mhd_flow.o
rm loader.txt
#
mv a.out ~/bin/$ARCH/mhd_flow
#
echo "Executable installed as ~/bin/$ARCH/mhd_flow"
