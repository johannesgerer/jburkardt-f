#!/bin/bash
#
gfortran -c -g random_data_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_data_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran random_data_prb.o -L$HOME/lib/$ARCH -lrandom_data
if [ $? -ne 0 ]; then
  echo "Errors linking and loading random_data_prb.o"
  exit
fi
rm random_data_prb.o
#
mv a.out random_data_prb
./random_data_prb > random_data_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running random_data_prb"
  exit
fi
rm random_data_prb
#
echo "Program output written to random_data_prb_output.txt"
