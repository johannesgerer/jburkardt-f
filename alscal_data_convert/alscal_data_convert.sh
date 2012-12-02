#!/bin/bash
#
gfortran -c -g alscal_data_convert.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling alscal_data_convert.f90"
  exit
fi
rm compiler.txt
#
gfortran alscal_data_convert.o
if [ $? -ne 0 ]; then
  echo "Errors loading alscal_data_convert.o"
  exit
fi
#
rm alscal_data_convert.o
mv a.out ~/bin/$ARCH/alscal_data_convert
#
echo "Executable installed as ~/bin/$ARCH/alscal_data_convert"
