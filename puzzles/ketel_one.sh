#!/bin/bash
#
gfortran -c -g ketel_one.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ketel_one.f90"
  exit
fi
rm compiler.txt
#
gfortran ketel_one.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ketel_one.o"
  exit
fi
rm ketel_one.o
#
mv a.out ketel_one
#
cp ../../datasets/words/wordlist.txt wordlist.txt
#
./ketel_one > ketel_one_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ketel_one"
  exit
fi
rm ketel_one
rm wordlist.txt
#
echo "Program output written to ketel_one_output.txt"
