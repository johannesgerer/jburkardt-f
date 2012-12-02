#!/bin/bash
#
gfortran -c -g dqed_prb3.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dqed_prb3.f90"
  exit
fi
rm compiler.txt
#
gfortran dqed_prb3.o -L$HOME/lib/$ARCH -ldqed
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dqed_prb3.o"
  exit
fi
rm dqed_prb3.o
#
mv a.out dqed_prb3
./dqed_prb3 > dqed_prb3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dqed_prb3"
  exit
fi
rm dqed_prb3
#
echo "Test program output written to dqed_prb3_output.txt."
