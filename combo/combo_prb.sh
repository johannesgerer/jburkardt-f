#!/bin/bash
#
gfortran -c -g combo_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling combo_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran combo_prb.o -L$HOME/lib/$ARCH -lcombo
if [ $? -ne 0 ]; then
  echo "Errors linking and loading combo_prb.o"
  exit
fi
rm combo_prb.o
#
mv a.out combo_prb
./combo_prb > combo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running combo_prb"
  exit
fi
rm combo_prb
#
echo "Test program output written to combo_prb_output.txt."
