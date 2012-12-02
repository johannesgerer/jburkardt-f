#!/bin/bash
#
gfortran -c -g machine_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machine_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran machine_prb.o -L$HOME/lib/$ARCH -lslatec
if [ $? -ne 0 ]; then
  echo "Errors linking and loading machine_prb.o"
  exit
fi
rm machine_prb.o
#
mv a.out machine_prb
./machine_prb > machine_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running machine_prb"
  exit
fi
rm machine_prb
#
echo "Program output written to machine_prb_output.txt"
