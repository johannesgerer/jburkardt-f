#!/bin/bash
#
gfortran -c -g xyz_demo.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xyz_demo.f90"
  exit
fi
rm compiler.txt
#
gfortran xyz_demo.o -L$HOME/lib/$ARCH -lplot3d_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xyz_demo.o"
  exit
fi
rm xyz_demo.o
#
mv a.out xyz_demo
./xyz_demo > xyz_demo_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running xyz_demo"
  exit
endif
rm xyz_demo
#
echo "Program output written to xyz_demo_output.txt"
