#!/bin/bash
#
gfortran -c -g pgma_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgma_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran pgma_io_prb.o -L$HOME/lib/$ARCH -lpgma_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pgma_io_prb.o"
  exit
fi
rm pgma_io_prb.o
#
mv a.out pgma_io_prb
./pgma_io_prb > pgma_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pgma_io_prb"
  exit
fi
rm pgma_io_prb
#
echo "Test program output written to pgma_io_prb_output.txt."
