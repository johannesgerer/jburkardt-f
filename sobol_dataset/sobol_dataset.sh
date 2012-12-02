#!/bin/bash
#
gfortran -c -g sobol_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sobol_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran sobol_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sobol_dataset.o"
  exit
fi
rm sobol_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/sobol_dataset
#
echo "Executable installed as ~/bin/$ARCH/sobol_dataset"
