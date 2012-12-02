#!/bin/bash
#
gfortran -c -g doomsday_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling doomsday_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran doomsday_prb.o -L$HOME/lib/$ARCH -ldoomsday
if [ $? -ne 0 ]; then
  echo "Errors linking and loading doomsday_prb.o"
  exit
fi
rm doomsday_prb.o
#
mv a.out doomsday_prb
./doomsday_prb > doomsday_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running doomsday_prb"
  exit
fi
rm doomsday_prb
#
echo "Program output written to doomsday_prb_output.txt"
