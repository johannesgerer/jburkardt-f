#!/bin/bash
#
gfortran -c -g select_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling select_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran select_prb.o -L$HOME/lib/$ARCH -lselect
if [ $? -ne 0 ]; then
  echo "Errors linking and loading select_prb.o"
  exit
fi
rm select_prb.o
#
mv a.out select_prb
./select_prb > select_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running select_prb"
  exit
fi
rm select_prb
#
echo "Test program output written to select_prb_output.txt."
