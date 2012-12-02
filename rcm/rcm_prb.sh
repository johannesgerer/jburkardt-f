#!/bin/bash
#
gfortran -c -g rcm_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rcm_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran rcm_prb.o -L$HOME/lib/$ARCH -lrcm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rcm_prb.o"
  exit
fi
rm rcm_prb.o
#
mv a.out rcm_prb
./rcm_prb > rcm_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rcm_prb"
  exit
fi
rm rcm_prb
#
echo "Test program output written to rcm_prb_output.txt."
