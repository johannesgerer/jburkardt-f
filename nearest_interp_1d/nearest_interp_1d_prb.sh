#!/bin/bash
#
gfortran -c -g nearest_interp_1d_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nearest_interp_1d_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran nearest_interp_1d_prb.o -L$HOME/lib/$ARCH -lnearest_interp_1d -ltest_interp -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nearest_interp_1d_prb.o"
  exit
fi
rm nearest_interp_1d_prb.o
#
mv a.out nearest_interp_1d_prb
./nearest_interp_1d_prb > nearest_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nearest_interp_1d_prb"
  exit
fi
rm nearest_interp_1d_prb
#
echo "Test program output written to nearest_interp_1d_prb_output.txt."
