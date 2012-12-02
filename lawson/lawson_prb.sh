#!/bin/bash
#
gfortran -c -g lawson_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lawson_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lawson_prb.o -L$HOME/lib/$ARCH -llawson
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lawson_prb.o"
  exit
fi
rm lawson_prb.o
#
mv a.out lawson_prb
./lawson_prb > lawson_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lawson_prb"
  exit
fi
rm lawson_prb
#
echo "Test program output written to lawson_prb_output.txt."
