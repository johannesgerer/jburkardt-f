#!/bin/bash
#
gfortran -c spacer_data_convert.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors spacer_data_convert.f90"
  exit
fi
rm compiler.txt
#
gfortran spacer_data_convert.o
if [ $? -ne 0 ]; then
  echo "Errors loading spacer_data_convert.o"
  exit
fi
#
rm spacer_data_convert.o
mv a.out ~/bin/$ARCH/spacer_data_convert
#
echo "Executable installed as ~/bin/$ARCH/spacer_data_convert"
