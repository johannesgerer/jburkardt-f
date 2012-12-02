#!/bin/bash
#
gfortran -c -g grafpack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grafpack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran grafpack_prb.o -L$HOME/lib/$ARCH -lgrafpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grafpack_prb.o"
  exit
fi
rm grafpack_prb.o
#
mv a.out grafpack_prb
./grafpack_prb > grafpack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running grafpack_prb"
  exit
fi
rm grafpack_prb
#
echo "Test program output written to grafpack_prb_output.txt."
