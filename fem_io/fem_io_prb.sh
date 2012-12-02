#!/bin/bash
#
gfortran -c -g fem_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran fem_io_prb.o -L$HOME/lib/$ARCH -lfem_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_io_prb.o"
  exit
fi
rm fem_io_prb.o
#
mv a.out fem_io_prb
./fem_io_prb > fem_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem_io_prb"
  exit
fi
rm fem_io_prb
#
echo "Test program output written to fem_io_prb_output.txt."
