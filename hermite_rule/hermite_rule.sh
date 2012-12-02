#!/bin/bash
#
gfortran -c -g hermite_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran hermite_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_rule.o"
  exit
fi
rm hermite_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/hermite_rule
#
echo "Executable installed as ~/bin/$ARCH/hermite_rule"
