#!/bin/bash
#
gfortran -c -g dqed_prb2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dqed_prb2.f90"
  exit
fi
rm compiler.txt
#
gfortran dqed_prb2.o -L$HOME/lib/$ARCH -ldqed
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dqed_prb2.o"
  exit
fi
rm dqed_prb2.o
#
mv a.out dqed_prb2
./dqed_prb2 > dqed_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dqed_prb2"
  exit
fi
rm dqed_prb2
#
echo "Test program output written to dqed_prb2_output.txt."
