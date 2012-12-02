#!/bin/bash
#
gfortran -c -g obj_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling obj_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran obj_io_prb.o -L$HOME/lib/$ARCH -lobj_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading obj_io_prb.o"
  exit
fi
rm obj_io_prb.o
#
mv a.out obj_io_prb
./obj_io_prb > obj_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running obj_io_prb"
  exit
fi
rm obj_io_prb
#
echo "Test program output written to obj_io_prb_output.txt."
