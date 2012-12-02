#!/bin/bash
#
gfortran -c -g eispack_prb1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling eispack_prb1.f90"
  exit
fi
rm compiler.txt
#
gfortran eispack_prb1.o -L$HOME/lib/$ARCH -leispack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading eispack_prb1.o"
  exit
fi
rm eispack_prb1.o
#
mv a.out eispack_prb1
./eispack_prb1 > eispack_prb1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running eispack_prb1"
  exit
fi
rm eispack_prb1
#
echo "Test program output written to eispack_prb1_output.txt."
