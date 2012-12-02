#!/bin/bash
#
gfortran -c -g dqed_prb4.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dqed_prb4.f90"
  exit
fi
rm compiler.txt
#
gfortran dqed_prb4.o -L$HOME/lib/$ARCH -ldqed
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dqed_prb4.o"
  exit
fi
rm dqed_prb4.o
#
mv a.out dqed_prb4
./dqed_prb4 > dqed_prb4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dqed_prb4"
  exit
fi
rm dqed_prb4
#
echo "Test program output written to dqed_prb4_output.txt."
