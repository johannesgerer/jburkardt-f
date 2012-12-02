#!/bin/bash
#
gfortran -c -g jacobi_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran jacobi_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading jacobi_rule.o"
  exit
fi
rm jacobi_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/jacobi_rule
#
echo "Executable installed as ~/bin/$ARCH/jacobi_rule"
