#!/bin/bash
#
gfortran -c airy_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling airy_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran airy_prb.o -L$HOME/lib/$ARCH -lslatec
if [ $? -ne 0 ]; then
  echo "Errors linking and loading airy_prb.o"
  exit
fi
rm airy_prb.o
#
mv a.out airy_prb
./airy_prb > airy_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running airy_prb"
  exit
fi
rm airy_prb
#
echo "Program output written to airy_prb_output.txt"
