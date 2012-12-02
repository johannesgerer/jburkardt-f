#!/bin/bash
#
gfortran -c -g nl2sol_prb2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nl2sol_prb2.f90"
  exit
fi
rm compiler.txt
#
gfortran nl2sol_prb2.o -L$HOME/lib/$ARCH -lnl2sol
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nl2sol_prb2.o"
  exit
fi
rm nl2sol_prb2.o
#
mv a.out nl2sol_prb2
./nl2sol_prb2 > nl2sol_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nl2sol_prb2"
  exit
fi
rm nl2sol_prb2
#
echo "Test program output written to nl2sol_prb2_output.txt."
