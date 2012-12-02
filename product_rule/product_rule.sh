#!/bin/bash
#
gfortran -c -g product_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling product_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran product_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading product_rule.o"
  exit
fi
rm product_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/product_rule
#
echo "Executable installed as ~/bin/$ARCH/product_rule"
