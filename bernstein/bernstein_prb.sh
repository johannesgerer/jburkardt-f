#!/bin/bash
#
gfortran -c -g bernstein_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bernstein_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran bernstein_prb.o -L$HOME/lib/$ARCH -lbernstein
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bernstein_prb.o"
  exit
fi
rm bernstein_prb.o
#
mv a.out bernstein_prb
./bernstein_prb > bernstein_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bernstein_prb"
  exit
fi
rm bernstein_prb
#
echo "Program output written to bernstein_prb_output.txt"
