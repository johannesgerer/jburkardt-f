#!/bin/bash
#
gfortran -c -g patterson_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling patterson_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran patterson_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading patterson_rule.o"
  exit
fi
rm patterson_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/patterson_rule
#
echo "Executable installed as ~/bin/$ARCH/patterson_rule"
