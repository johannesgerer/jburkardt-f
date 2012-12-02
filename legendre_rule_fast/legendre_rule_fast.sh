#!/bin/bash
#
gfortran -c -g legendre_rule_fast.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_rule_fast.f90"
  exit
fi
rm compiler.txt
#
gfortran legendre_rule_fast.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading legendre_rule_fast.o"
  exit
fi
rm legendre_rule_fast.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/legendre_rule_fast
#
echo "Executable installed as ~/bin/$ARCH/legendre_rule_fast"
