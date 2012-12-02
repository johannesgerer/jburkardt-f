#!/bin/bash
#
gfortran -c -g box_behnken_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling box_behnken_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran box_behnken_prb.o -L$HOME/lib/$ARCH -lbox_behnken
if [ $? -ne 0 ]; then
  echo "Errors linking and loading box_behnken_prb.o"
  exit
fi
rm box_behnken_prb.o
#
mv a.out box_behnken_prb
./box_behnken_prb > box_behnken_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running box_behnken_prb"
  exit
fi
rm box_behnken_prb
#
echo "Test program output written to box_behnken_prb_output.txt."
