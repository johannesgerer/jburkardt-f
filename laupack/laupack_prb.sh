#!/bin/bash
#
gfortran -c -g laupack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laupack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran laupack_prb.o -L$HOME/lib/$ARCH -llaupack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laupack_prb.o"
  exit
fi
rm laupack_prb.o
#
mv a.out laupack_prb
./laupack_prb > laupack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running laupack_prb"
  exit
fi
rm laupack_prb
#
echo "Test program output written to laupack_prb_output.txt."
