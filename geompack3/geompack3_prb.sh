#!/bin/bash
#
gfortran -c -g geompack3_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geompack3_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran geompack3_prb.o -L$HOME/lib/$ARCH -lgeompack3
if [ $? -ne 0 ]; then
  echo "Errors linking and loading geompack3_prb.o"
  exit
fi
rm geompack3_prb.o
#
mv a.out geompack3_prb
./geompack3_prb > geompack3_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running geompack3_prb"
  exit
fi
rm geompack3_prb
#
echo "Test program output written to geompack3_prb_output.txt."
