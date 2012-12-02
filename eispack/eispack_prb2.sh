#!/bin/bash
#
gfortran -c -g eispack_prb2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling eispack_prb2.f90"
  exit
fi
rm compiler.txt
#
gfortran eispack_prb2.o -L$HOME/lib/$ARCH -leispack -ltest_eigen
if [ $? -ne 0 ]; then
  echo "Errors linking and loading eispack_prb2.o"
  exit
fi
rm eispack_prb2.o
#
mv a.out eispack_prb2
./eispack_prb2 > eispack_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running eispack_prb2"
  exit
fi
rm eispack_prb2
#
echo "Test program output written to eispack_prb2_output.txt."
