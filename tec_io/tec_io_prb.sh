#!/bin/bash
#
gfortran -c -g tec_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tec_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran tec_io_prb.o -L$HOME/lib/$ARCH -ltec_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tec_io_prb.o"
  exit
fi
rm tec_io_prb.o
#
mv a.out tec_io_prb
./tec_io_prb > tec_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tec_io_prb"
  exit
fi
rm tec_io_prb
#
echo "Test program output written to tec_io_prb_output.txt."
