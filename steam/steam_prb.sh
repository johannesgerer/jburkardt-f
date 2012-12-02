#!/bin/bash
#
gfortran -c -g steam_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling steam_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran steam_prb.o -L$HOME/lib/$ARCH -lsteam
if [ $? -ne 0 ]; then
  echo "Errors linking and loading steam_prb.o"
  exit
fi
rm steam_prb.o
#
mv a.out steam_prb
./steam_prb > steam_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running steam_prb"
  exit
fi
rm steam_prb
#
echo "Test program output written to steam_prb_output.txt."
