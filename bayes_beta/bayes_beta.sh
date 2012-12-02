#!/bin/bash
#
gfortran -c bayes_beta.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bayes_beta.f90"
  exit
fi
rm compiler.txt
#
gfortran bayes_beta.o
if [ $? -ne 0 ]; then
  echo "Errors loading bayes_beta.o"
  exit
fi
rm bayes_beta.o
#
mv a.out ~/bin/$ARCH/bayes_beta
#
echo "Executable installed as ~/bin/$ARCH/bayes_beta"
