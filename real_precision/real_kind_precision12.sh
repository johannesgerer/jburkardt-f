#!/bin/bash
#
gfortran -c -g real_kind_precision12.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling real_kind_precision12.f90"
  exit
fi
rm compiler.txt
#
gfortran real_kind_precision12.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading real_kind_precision12.o"
  exit
fi
rm real_kind_precision12.o
#
mv a.out real_kind_precision12
./real_kind_precision12 > real_kind_precision12_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running real_kind_precision12"
  exit
fi
rm real_kind_precision12
#
echo "The real_kind_precision12 test problem has been executed."
