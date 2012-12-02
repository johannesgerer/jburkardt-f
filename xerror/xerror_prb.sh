#!/bin/bash
#
gfortran -c -g xerror_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xerror_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran xerror_prb.o -L$HOME/lib/$ARCH -lxerror
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xerror_prb.o"
  exit
fi
rm xerror_prb.o
#
mv a.out xerror_prb
./xerror_prb > xerror_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running xerror_prb"
  exit
fi
rm xerror_prb
#
echo "Test program output written to xerror_prb_output.txt."
