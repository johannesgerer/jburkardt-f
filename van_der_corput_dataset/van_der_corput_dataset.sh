#!/bin/bash
#
gfortran -c -g van_der_corput_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling van_der_corput_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran van_der_corput_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading van_der_corput_dataset.o"
  exit
fi
rm van_der_corput_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/van_der_corput_dataset
#
echo "Executable installed as ~/bin/$ARCH/van_der_corput_dataset"
