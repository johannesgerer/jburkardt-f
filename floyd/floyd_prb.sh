#!/bin/bash
#
gfortran -c -g floyd_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling floyd_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran floyd_prb.o -L$HOME/lib/$ARCH -lfloyd
if [ $? -ne 0 ]; then
  echo "Errors linking and loading floyd_prb.o"
  exit
fi
rm floyd_prb.o
#
mv a.out floyd_prb
./floyd_prb > floyd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running floyd_prb"
  exit
fi
rm floyd_prb
#
echo "Test program output written to floyd_prb_output.txt."
