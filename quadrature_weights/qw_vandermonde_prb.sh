#!/bin/bash
#
gfortran -c -g qw_vandermonde_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qw_vandermonde_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran qw_vandermonde_prb.o -L$HOME/lib/$ARCH -lqw_vandermonde
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qw_vandermonde_prb.o"
  exit
fi
rm qw_vandermonde_prb.o
#
mv a.out qw_vandermonde_prb
./qw_vandermonde_prb > qw_vandermonde_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qw_vandermonde_prb"
  exit
fi
rm qw_vandermonde_prb
#
echo "Program output written to qw_vandermonde_prb_output.txt"
