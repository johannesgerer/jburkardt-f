#!/bin/bash
#
gfortran -c -g linpack_s_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_s_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran linpack_s_prb.o -L$HOME/lib/$ARCH -llinpack_s -lblas1_s
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_s_prb.o"
  exit
fi
rm linpack_s_prb.o
#
mv a.out linpack_s_prb
./linpack_s_prb > linpack_s_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linpack_s_prb"
  exit
fi
rm linpack_s_prb
#
echo "Test program output written to linpack_s_prb_output.txt."
