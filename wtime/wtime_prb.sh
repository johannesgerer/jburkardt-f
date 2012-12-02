#!/bin/bash
#
gfortran -c -g wtime_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran wtime_prb.o -L$HOME/lib/$ARCH -lwtime
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wtime_prb.o"
  exit
fi
rm wtime_prb.o
#
mv a.out wtime_prb
./wtime_prb > wtime_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wtime_prb"
  exit
fi
rm wtime_prb
#
echo "Test program output written to wtime_prb_output.txt."
