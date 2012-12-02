#!/bin/bash
#
gfortran -c -g gegenbauer_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gegenbauer_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran gegenbauer_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gegenbauer_rule.o"
  exit
fi
rm gegenbauer_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/gegenbauer_rule
#
echo "Executable installed as ~/bin/$ARCH/gegenbauer_rule"
