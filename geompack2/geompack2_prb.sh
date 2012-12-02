#!/bin/bash
#
gfortran -c -g geompack2_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geompack2_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran geompack2_prb.o -L$HOME/lib/$ARCH -lgeompack2
if [ $? -ne 0 ]; then
  echo "Errors linking and loading geompack2_prb.o"
  exit
fi
rm geompack2_prb.o
#
mv a.out geompack2_prb
./geompack2_prb > geompack2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running geompack2_prb"
  exit
fi
rm geompack2_prb
#
echo "Test program output written to geompack2_prb_output.txt."
