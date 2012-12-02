#!/bin/bash
#
gfortran -c -g vector_max.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vector_max.f90"
  exit
fi
rm compiler.txt
#
gfortran vector_max.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vector_max.o"
  exit
fi
rm vector_max.o
#
mv a.out vector_max
./vector_max > vector_max_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running vector_max"
  exit
fi
rm vector_max
#
echo "Program output written to vector_max_output.txt"
