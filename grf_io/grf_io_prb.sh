#!/bin/bash
#
gfortran -c -g grf_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grf_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran grf_io_prb.o -L$HOME/lib/$ARCH -lgrf_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grf_io_prb.o"
  exit
fi
rm grf_io_prb.o
#
mv a.out grf_io_prb
./grf_io_prb > grf_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running grf_io_prb"
  exit
fi
rm grf_io_prb
#
echo "Test program output written to grf_io_prb_output.txt."
