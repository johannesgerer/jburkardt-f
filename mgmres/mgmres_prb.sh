#!/bin/bash
#
gfortran -c -g mgmres_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mgmres_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran mgmres_prb.o -L$HOME/lib/$ARCH -lmgmres
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mgmres_prb.o"
  exit
fi
rm mgmres_prb.o
#
mv a.out mgmres_prb
./mgmres_prb > mgmres_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mgmres_prb"
  exit
fi
rm mgmres_prb
#
echo "Test program output written to mgmres_prb_output.txt."
