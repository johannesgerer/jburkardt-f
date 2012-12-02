#!/bin/bash
#
gfortran -c -g real_kind_precision18.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling real_kind_precision18.f90"
  exit
fi
rm compiler.txt
#
gfortran real_kind_precision18.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading real_kind_precision18.o"
  exit
fi
rm real_kind_precision18.o
#
mv a.out real_kind_precision18
./real_kind_precision18 > real_kind_precision18_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running real_kind_precision18"
  exit
fi
rm real_kind_precision18
#
echo "The real_kind_precision18 test problem has been executed."
