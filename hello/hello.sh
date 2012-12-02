#!/bin/bash
#
gfortran -c -g hello.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hello.f90"
  exit
fi
rm compiler.txt
#
gfortran hello.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hello.o"
  exit
fi
rm hello.o
#
mv a.out hello
./hello > hello_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hello"
  exit
fi
rm hello
#
echo "Program output written to hello_output.txt"
