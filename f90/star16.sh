#!/bin/bash
#
g95 -c -g star16.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling star16.f90"
  exit
fi
rm compiler.txt
#
g95 star16.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading star16.o"
  exit
fi
rm star16.o
#
mv a.out star16
./star16 > star16_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running star16"
  exit
fi
rm star16
#
echo "The star16 test problem has been executed."
