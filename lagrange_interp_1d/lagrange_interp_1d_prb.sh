#!/bin/bash
#
gfortran -c -g lagrange_interp_1d_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_interp_1d_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lagrange_interp_1d_prb.o -L$HOME/lib/$ARCH -llagrange_interp_1d -ltest_interp_1d -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lagrange_interp_1d_prb.o"
  exit
fi
rm lagrange_interp_1d_prb.o
#
mv a.out lagrange_interp_1d_prb
./lagrange_interp_1d_prb > lagrange_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lagrange_interp_1d_prb"
  exit
fi
rm lagrange_interp_1d_prb
#
echo "Test program output written to lagrange_interp_1d_prb_output.txt."
