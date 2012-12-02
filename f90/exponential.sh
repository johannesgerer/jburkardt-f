#!/bin/bash
#
gfortran -c -g exponential.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling exponential.f90"
  exit
fi
rm compiler.txt
#
gfortran exponential.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading exponential.o"
  exit
fi
rm exponential.o
#
mv a.out exponential
./exponential > exponential_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running exponential"
  exit
fi
rm exponential
#
echo "The exponential test problem has been executed."
