#!/bin/bash
#
gfortran -c -g dlap_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dlap_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran dlap_io_prb.o -L$HOME/lib/$ARCH -ldlap_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dlap_io_prb.o"
  exit
fi
rm dlap_io_prb.o
#
mv a.out dlap_io_prb
./dlap_io_prb > dlap_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dlap_io_prb"
  exit
fi
rm dlap_io_prb
#
echo "Test program output written to dlap_io_prb_output.txt."
