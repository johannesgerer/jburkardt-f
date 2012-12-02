#!/bin/bash
#
gfortran -c -g poly_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poly_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran poly_io_prb.o -L$HOME/lib/$ARCH -lpoly_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading poly_io_prb.o"
  exit
fi
rm poly_io_prb.o
#
mv a.out poly_io_prb
./poly_io_prb > poly_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poly_io_prb"
  exit
fi
rm poly_io_prb
#
echo "Test program output written to poly_io_prb_output.txt."
