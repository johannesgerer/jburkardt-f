#!/bin/bash
#
gfortran -c -g unigram.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling unigram.f90"
  exit
fi
rm compiler.txt
#
gfortran unigram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading unigram.o"
  exit
fi
rm unigram.o
#
mv a.out unigram
#
cp ../../datasets/words/wordlist.txt wordlist.txt
#
./unigram > unigram_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running unigram"
  exit
fi
rm unigram
rm wordlist.txt
#
echo "Program output written to unigram_output.txt"
