#!/bin/bash
#
gfortran -c -g grid_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grid_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran grid_prb.o -L$HOME/lib/$ARCH -lgrid
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grid_prb.o"
  exit
fi
rm grid_prb.o
#
mv a.out grid_prb
./grid_prb > grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running grid_prb"
  exit
fi
rm grid_prb
#
echo "Test program output written to grid_prb_output.txt."
