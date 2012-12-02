#!/bin/bash
#
gfortran -c -g pitcon7_prb5.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pitcon7_prb5.f90"
  exit
fi
rm compiler.txt
#
gfortran pitcon7_prb5.o -L$HOME/lib/$ARCH -lpitcon7
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pitcon7_prb5.o"
  exit
fi
rm pitcon7_prb5.o
#
mv a.out pitcon7_prb5
./pitcon7_prb5 > pitcon7_prb5_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pitcon7_prb5"
  exit
fi
rm pitcon7_prb5
#
echo "Test program output written to pitcon7_prb5_output.txt."
