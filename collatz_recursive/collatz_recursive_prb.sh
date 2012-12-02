#!/bin/bash
#
gfortran -c -g collatz_recursive_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling collatz_recursive_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran collatz_recursive_prb.o -L$HOME/lib/$ARCH -lcollatz_recursive
if [ $? -ne 0 ]; then
  echo "Errors linking and loading collatz_recursive_prb.o"
  exit
fi
rm collatz_recursive_prb.o
#
mv a.out collatz_recursive_prb
./collatz_recursive_prb > collatz_recursive_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running collatz_recursive_prb"
  exit
fi
rm collatz_recursive_prb
#
echo "Test program output written to collatz_recursive_prb_output.txt."
