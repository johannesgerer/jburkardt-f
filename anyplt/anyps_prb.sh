#!/bin/bash
#
gfortran -c -g anyplt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran anyplt_prb.o -L$HOME/lib/$ARCH -lanyps -lps_write
if [ $? -ne 0 ]; then
  echo "Errors linking and loading anyplt_prb.o"
  exit
fi
rm anyplt_prb.o
#
mv a.out anyps_prb
./anyps_prb > anyps_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running anyps_prb"
fi
rm anyps_prb
#
echo "Program output written to anyps_prb_output.txt"
