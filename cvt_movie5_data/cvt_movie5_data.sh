#!/bin/bash
#
gfortran -c -g cvt_movie5_data.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_movie5_data.csh"
  exit
fi
rm compiler.txt
#
gfortran cvt_movie5_data.o -L$HOME/lib/$ARCH -ltest_triangulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_movie5_data.o"
  exit
fi
rm cvt_movie5_data.o
#
mv a.out cvt_movie5_data
./cvt_movie5_data > cvt_movie5_data_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_movie5_data"
  exit
fi
rm cvt_movie5_data
#
echo "Program output written to cvt_movie5_data_output.txt"
