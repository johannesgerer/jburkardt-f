#!/bin/bash
#
gfortran -c -g niederreiter2_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter2_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran niederreiter2_prb.o -L$HOME/lib/$ARCH -lniederreiter2
if [ $? -ne 0 ]; then
  echo "Errors linking and loading niederreiter2_prb.o"
  exit
fi
rm niederreiter2_prb.o
#
mv a.out niederreiter2_prb
./niederreiter2_prb > niederreiter2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running niederreiter2_prb"
  exit
fi
rm niederreiter2_prb
#
echo "Test program output written to niederreiter2_prb_output.txt."
