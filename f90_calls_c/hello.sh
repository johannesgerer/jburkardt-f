#!/bin/bash
#
#  Compile the main FORTRAN90 program.
#
gfortran -c -g hello_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hello_prb.f90"
  exit
fi
rm compiler.txt
#
#  Compile the C function.
#
gcc -c -g hello.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hello_prb.c"
  exit
fi
rm compiler.txt
#
#  Link the Main program and the function.
#
gfortran hello_prb.o hello.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hello_prb.o + hello.o"
  exit
fi
rm hello_prb.o
rm hello.o
#
#  Run the program.
#
mv a.out hello
./hello > hello_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hello"
  exit
fi
rm hello
#
echo "Test program output written to hello_output.txt."
