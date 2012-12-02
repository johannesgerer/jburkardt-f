#!/bin/bash
#
gfortran -c -g index_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling index_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran index_prb.o -L$HOME/lib/$ARCH -lindex
if [ $? -ne 0 ]; then
  echo "Errors linking and loading index_prb.o"
  exit
fi
rm index_prb.o
#
mv a.out index_prb
./index_prb > index_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running index_prb"
  exit
fi
rm index_prb
#
echo "Test program output written to index_prb_output.txt."
