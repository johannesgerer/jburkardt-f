#!/bin/bash
#
gfortran -c -g laguerre_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran laguerre_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laguerre_rule.o"
  exit
fi
rm laguerre_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/laguerre_rule
#
echo "Executable installed as ~/bin/$ARCH/laguerre_rule"
