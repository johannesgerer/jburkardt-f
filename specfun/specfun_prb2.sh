#!/bin/bash
#
gfortran -c -g specfun_prb2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling specfun_prb2.f90"
  exit
fi
rm compiler.txt
#
gfortran specfun_prb2.o -L$HOME/lib/$ARCH -lspecfun
if [ $? -ne 0 ]; then
  echo "Errors linking and loading specfun_prb2.o"
  exit
fi
rm specfun_prb2.o
#
mv a.out specfun_prb2
./specfun_prb2 > specfun_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running specfun_prb2"
  exit
fi
rm specfun_prb2
#
echo "Test program output written to specfun_prb2_output.txt."
