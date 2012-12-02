#!/bin/bash
#
gfortran -c bayes_weight.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bayes_weight.f90"
  exit
fi
rm compiler.txt
#
gfortran bayes_weight.o
if [ $? -ne 0 ]; then
  echo "Errors loading bayes_weight.o"
  exit
fi
rm bayes_weight.o
#
mv a.out ~/bin/$ARCH/bayes_weight
#
echo "Executable installed as ~/bin/$ARCH/bayes_weight"
