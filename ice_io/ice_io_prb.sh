#!/bin/bash
#
g95 -c -g ice_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ice_io_prb.f90"
  exit
fi
rm compiler.txt
#
g95 ice_io_prb.o -L$HOME/lib/$ARCH -lice_io -L/usr/local/lib -lnetcdf
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ice_io_prb.o"
  exit
fi
rm ice_io_prb.o
#
mv a.out ice_io_prb
./ice_io_prb > ice_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ice_io_prb"
  exit
fi
rm ice_io_prb
#
echo "Test program output written to ice_io_prb_output.txt."
