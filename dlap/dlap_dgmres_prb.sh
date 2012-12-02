#!/bin/bash
#
gfortran -c -g dlap_dgmres_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dlap_dgmres_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran dlap_dgmres_prb.o -L$HOME/lib/$ARCH -ldlap
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dlap_dgmres_prb.o"
  exit
fi
rm dlap_dgmres_prb.o
#
mv a.out dlap_dgmres_prb
./dlap_dgmres_prb > dlap_dgmres_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dlap_dgmres_prb"
  exit
fi
rm dlap_dgmres_prb
#
echo "Test program output written to dlap_dgmres_prb_output.txt."
