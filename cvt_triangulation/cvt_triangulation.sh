#!/bin/bash
#
gfortran -c -g cvt_triangulation.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_triangulation.csh"
  exit
fi
rm compiler.txt
#
gfortran cvt_triangulation.o -L$HOME/lib/$ARCH -ltest_triangulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_triangulation.o"
  exit
fi
rm cvt_triangulation.o
#
mv a.out cvt_triangulation
./cvt_triangulation > cvt_triangulation_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_triangulation"
  exit
fi
rm cvt_triangulation
#
echo "Program output written to cvt_triangulation_output.txt"
