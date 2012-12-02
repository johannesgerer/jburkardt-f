#!/bin/bash
#
gfortran -c -g mxm.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxm.f90"
  exit
fi
rm compiler.txt
#
gfortran mxm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxm.o"
  exit
fi
rm mxm.o
#
mv a.out mxm
./mxm > mxm_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mxm"
  exit
fi
rm mxm
#
echo "The mxm test problem has been executed."
