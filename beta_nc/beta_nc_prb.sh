#!/bin/bash
#
gfortran -c -g beta_nc_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling beta_nc_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran beta_nc_prb.o -L$HOME/lib/$ARCH -lbeta_nc
if [ $? -ne 0 ]; then
  echo "Errors linking and loading beta_nc_prb.o"
  exit
fi
rm beta_nc_prb.o
#
mv a.out beta_nc_prb
./beta_nc_prb > beta_nc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running beta_nc_prb"
  exit
fi
rm beta_nc_prb
#
echo "Test program output written to beta_nc_prb_output.txt."
