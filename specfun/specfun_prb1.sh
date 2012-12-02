#!/bin/bash
#
gfortran -c -g specfun_prb1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling specfun_prb1.f90"
  exit
fi
rm compiler.txt
#
gfortran specfun_prb1.o -L$HOME/lib/$ARCH -lspecfun
if [ $? -ne 0 ]; then
  echo "Errors linking and loading specfun_prb1.o"
  exit
fi
rm specfun_prb1.o
#
mv a.out specfun_prb1
./specfun_prb1 > specfun_prb1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running specfun_prb1"
  exit
fi
rm specfun_prb1
#
echo "Test program output written to specfun_prb1_output.txt."
