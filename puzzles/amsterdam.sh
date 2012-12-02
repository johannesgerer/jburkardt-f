#!/bin/bash
#
gfortran -c -g amsterdam.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling amsterdam.f90"
  exit
fi
rm compiler.txt
#
gfortran amsterdam.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading amsterdam.o"
  exit
fi
rm amsterdam.o
#
mv a.out amsterdam
#
cp ../../datasets/words/wordlist.txt wordlist.txt
#
./amsterdam > amsterdam_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running amsterdam"
  exit
fi
rm amsterdam
rm wordlist.txt
#
echo "Program output written to amsterdam_output.txt"
