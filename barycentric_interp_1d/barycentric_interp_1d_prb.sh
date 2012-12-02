#!/bin/bash
#
gfortran -c -g barycentric_interp_1d_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling barycentric_interp_1d_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran barycentric_interp_1d_prb.o -L$HOME/lib/$ARCH -lbarycentric_interp_1d \
  -ltest_interp_1d -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading barycentric_interp_1d_prb.o"
  exit
fi
rm barycentric_interp_1d_prb.o
#
mv a.out barycentric_interp_1d_prb
./barycentric_interp_1d_prb > barycentric_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running barycentric_interp_1d_prb"
  exit
fi
rm barycentric_interp_1d_prb
#
echo "Test program output written to barycentric_interp_1d_prb_output.txt."
