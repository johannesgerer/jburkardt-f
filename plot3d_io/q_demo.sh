#!/bin/bash
#
gfortran -c -g q_demo.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling q_demo.f90"
  exit
fi
rm compiler.txt
#
gfortran q_demo.o -L$HOME/lib/$ARCH -lplot3d_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading q_demo.o"
  exit
fi
rm q_demo.o
#
mv a.out q_demo
./q_demo > q_demo_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running q_demo"
  exit
fi
rm q_demo
#
echo "Program output written to q_demo_output.txt"
