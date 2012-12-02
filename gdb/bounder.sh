#!/bin/bash
#
gfortran -c -g bounder.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bounder.f90"
  exit
fi
rm compiler.txt
#
gfortran bounder.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bounder.o"
  exit
fi
rm bounder.o
#
mv a.out bounder
./bounder > bounder_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bounder"
  exit
fi
rm bounder
#
echo "Program output written to bounder_output.txt"
