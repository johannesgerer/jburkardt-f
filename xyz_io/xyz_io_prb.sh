#!/bin/bash
#
gfortran -c -g xyz_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xyz_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran xyz_io_prb.o -L$HOME/lib/$ARCH -lxyz_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xyz_io_prb.o"
  exit
fi
rm xyz_io_prb.o
#
mv a.out xyz_io_prb
./xyz_io_prb > xyz_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running xyz_io_prb"
  exit
fi
rm xyz_io_prb
#
echo "Test program output written to xyz_io_prb_output.txt."
