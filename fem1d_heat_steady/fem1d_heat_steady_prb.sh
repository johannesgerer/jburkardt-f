#!/bin/bash
#
gfortran -c -g fem1d_heat_steady_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_heat_steady_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran fem1d_heat_steady_prb.o -L$HOME/lib/$ARCH -lfem1d_heat_steady
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_heat_steady_prb.o"
  exit
fi
rm fem1d_heat_steady_prb.o
#
mv a.out fem1d_heat_steady_prb
./fem1d_heat_steady_prb > fem1d_heat_steady_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem1d_heat_steady_prb"
  exit
fi
rm fem1d_heat_steady_prb
#
echo "Program output written to fem1d_heat_steady_prb_output.txt"
