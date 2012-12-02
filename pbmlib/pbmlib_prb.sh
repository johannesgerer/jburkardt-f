#!/bin/bash
#
gfortran -c -g pbmlib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbmlib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran pbmlib_prb.o -L$HOME/lib/$ARCH -lpbmlib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pbmlib_prb.o"
  exit
fi
rm pbmlib_prb.o
#
mv a.out pbmlib_prb
./pbmlib_prb > pbmlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pbmlib_prb"
  exit
fi
rm pbmlib_prb
#
echo "Test program output written to pbmlib_prb_output.txt."
