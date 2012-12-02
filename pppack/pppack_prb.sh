#!/bin/bash
#
gfortran -c -g pppack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pppack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran pppack_prb.o -L$HOME/lib/$ARCH -lpppack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pppack_prb.o"
  exit
fi
rm pppack_prb.o
#
mv a.out pppack_prb
./pppack_prb > pppack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pppack_prb"
  exit
fi
rm pppack_prb
#
echo "Test program output written to pppack_prb_output.txt."
