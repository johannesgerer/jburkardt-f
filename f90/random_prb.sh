#!/bin/bash
#
gfortran -c -g random_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran random_prb.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading random_prb.o"
  exit
fi
rm random_prb.o
#
mv a.out random_prb
./random_prb > random_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running random_prb"
  exit
fi
rm random_prb
#
echo "The random_prb test problem has been executed."
