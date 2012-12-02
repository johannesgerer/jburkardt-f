#!/bin/bash
#
gfortran -c -g real_kind_default.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling real_kind_default.f90"
  exit
fi
rm compiler.txt
#
gfortran real_kind_default.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading real_kind_default.o"
  exit
fi
rm real_kind_default.o
#
mv a.out real_kind_default
./real_kind_default > real_kind_default_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running real_kind_default"
  exit
fi
rm real_kind_default
#
echo "The real_kind_default test problem has been executed."
