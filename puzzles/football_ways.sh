#!/bin/bash
#
gfortran -c -g football_ways.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling football_ways.f90"
  exit
fi
rm compiler.txt
#
gfortran football_ways.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading football_ways.o"
  exit
fi
rm football_ways.o
#
mv a.out football_ways
#
./football_ways > football_ways_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running football_ways"
  exit
fi
rm football_ways
#
echo "Program output written to football_ways_output.txt"
