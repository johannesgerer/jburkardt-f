#!/bin/bash
#
gfortran -c -g csv_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling csv_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran csv_io_prb.o -L$HOME/lib/$ARCH -lcsv_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading csv_io_prb.o"
  exit
fi
rm csv_io_prb.o
#
mv a.out csv_io_prb
./csv_io_prb > csv_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running csv_io_prb"
  exit
fi
rm csv_io_prb
#
echo "Test program output written to csv_io_prb_output.txt."
