#!/bin/bash
#
gfortran -c -g clenshaw_curtis_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling clenshaw_curtis_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran clenshaw_curtis_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clenshaw_curtis_rule.o"
  exit
fi
rm clenshaw_curtis_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/clenshaw_curtis_rule
#
echo "Executable installed as ~/bin/$ARCH/clenshaw_curtis_rule"
