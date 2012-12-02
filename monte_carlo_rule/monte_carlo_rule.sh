#!/bin/bash
#
gfortran -c -g monte_carlo_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling monte_carlo_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran monte_carlo_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading monte_carlo_rule.o"
  exit
fi
rm monte_carlo_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/monte_carlo_rule
#
echo "Executable installed as ~/bin/$ARCH/monte_carlo_rule"
