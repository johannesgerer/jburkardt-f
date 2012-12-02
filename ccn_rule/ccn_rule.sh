#!/bin/bash
#
gfortran -c -g ccn_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ccn_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran ccn_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ccn_rule.o"
  exit
fi
rm ccn_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/ccn_rule
#
echo "Executable installed as ~/bin/$ARCH/ccn_rule"
