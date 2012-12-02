#!/bin/bash
#
gfortran -c -g rejoin_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rejoin_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran rejoin_prb.o -L$HOME/lib/$ARCH -lrejoin
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rejoin_prb.o"
  exit
fi
rm rejoin_prb.o
#
mv a.out rejoin_prb
./rejoin_prb > rejoin_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rejoin_prb"
  exit
fi
rm rejoin_prb
#
echo "Test program output written to rejoin_prb_output.txt."
