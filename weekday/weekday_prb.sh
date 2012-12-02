#!/bin/bash
#
gfortran -c -g weekday_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling weekday_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran weekday_prb.o -L$HOME/lib/$ARCH -lweekday
if [ $? -ne 0 ]; then
  echo "Errors linking and loading weekday_prb.o"
  exit
fi
rm weekday_prb.o
#
mv a.out weekday_prb
./weekday_prb > weekday_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running weekday_prb"
  exit
fi
rm weekday_prb
#
echo "Test program output written to weekday_prb_output.txt."
