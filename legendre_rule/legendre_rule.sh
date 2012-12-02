#!/bin/bash
#
gfortran -c -g legendre_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran legendre_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading legendre_rule.o"
  exit
fi
rm legendre_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/legendre_rule
#
echo "Executable installed as ~/bin/$ARCH/legendre_rule"
