#!/bin/bash
#
ar x $HOME/lib/$ARCH/libhb_io.a hb_file_module.mod
#
gfortran -c -g hb_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hb_io_prb.f90"
  exit
fi
rm compiler.txt
rm hb_file_module.mod
#
gfortran hb_io_prb.o -L$HOME/lib/$ARCH -lhb_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hb_io_prb.o"
  exit
fi
rm hb_io_prb.o
#
mv a.out hb_io_prb
./hb_io_prb > hb_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hb_io_prb"
  exit
fi
rm hb_io_prb
#
echo "Test program output written to hb_io_prb_output.txt."
