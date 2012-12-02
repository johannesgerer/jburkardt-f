#!/bin/bash
#
gfortran -c -g md.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md.f90"
  exit
fi
rm compiler.txt
#
gfortran md.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md.o"
  exit
fi
rm md.o
#
mv a.out md
./md < md_input.txt > md_output.txt
rm md
#
echo "Output written to md_output.txt"
