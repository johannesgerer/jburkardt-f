#!/bin/bash
#
gfortran -c -g gen_laguerre_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gen_laguerre_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran gen_laguerre_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gen_laguerre_rule.o"
  exit
fi
rm gen_laguerre_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/gen_laguerre_rule
#
echo "Executable installed as ~/bin/$ARCH/gen_laguerre_rule"
