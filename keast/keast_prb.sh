#!/bin/bash
#
gfortran -c -g keast_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling keast_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran keast_prb.o -L$HOME/lib/$ARCH -lkeast
if [ $? -ne 0 ]; then
  echo "Errors linking and loading keast_prb.o"
  exit
fi
rm keast_prb.o
#
mv a.out keast_prb
./keast_prb > keast_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running keast_prb"
  exit
fi
rm keast_prb
#
echo "Test program output written to keast_prb_output.txt."
