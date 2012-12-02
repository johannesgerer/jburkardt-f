#!/bin/bash
#
gfortran -c -g power_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling power_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran power_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading power_rule.o"
  exit
fi
rm power_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/power_rule
#
echo "Executable installed as ~/bin/$ARCH/power_rule"
