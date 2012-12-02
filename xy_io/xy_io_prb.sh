#!/bin/bash
#
gfortran -c -g xy_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xy_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran xy_io_prb.o -L$HOME/lib/$ARCH -lxy_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xy_io_prb.o"
  exit
fi
rm xy_io_prb.o
#
mv a.out xy_io_prb
./xy_io_prb > xy_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running xy_io_prb"
  exit
fi
rm xy_io_prb
#
echo "Test program output written to xy_io_prb_output.txt."
