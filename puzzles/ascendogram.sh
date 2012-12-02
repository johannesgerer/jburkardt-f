#!/bin/bash
#
gfortran -c -g ascendogram.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ascendogram.f90"
  exit
fi
rm compiler.txt
#
gfortran ascendogram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ascendogram.o"
  exit
fi
rm ascendogram.o
#
mv a.out ascendogram
#
cp ../../datasets/words/wordlist.txt wordlist.txt
#
./ascendogram > ascendogram_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ascendogram"
  exit
fi
rm ascendogram
rm wordlist.txt
#
echo "Program output written to ascendogram_output.txt"
