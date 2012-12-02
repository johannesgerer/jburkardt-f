#!/bin/bash
#
gfortran -c -g niederreiter2_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter2_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran niederreiter2_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading niederreiter2_dataset.o"
  exit
fi
rm niederreiter2_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/niederreiter2_dataset
#
echo "Executable installed as ~/bin/$ARCH/niederreiter2_dataset"
