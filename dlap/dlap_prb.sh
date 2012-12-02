#!/bin/bash
#
gfortran -c -g dlap_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dlap_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran dlap_prb.o -L$HOME/lib/$ARCH -ldlap
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dlap_prb.o"
  exit
fi
rm dlap_prb.o
#
mv a.out dlap_prb
./dlap_prb > dlap_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dlap_prb"
  exit
fi
rm dlap_prb
#
echo "Test program output written to dlap_prb_output.txt."
