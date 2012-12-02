#!/bin/bash
#
gfortran -c -g tanh_quad_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tanh_quad_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran tanh_quad_prb.o -L$HOME/lib/$ARCH -ltanh_quad
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tanh_quad_prb.o"
  exit
fi
rm tanh_quad_prb.o
#
mv a.out tanh_quad_prb
./tanh_quad_prb > tanh_quad_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tanh_quad_prb"
  exit
fi
rm tanh_quad_prb
#
echo "Test program output written to tanh_quad_prb_output.txt."
