#!/bin/bash
#
gfortran -c -g basis_compare.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling basis_compare.f90"
  exit
fi
rm compiler.txt
#
gfortran basis_compare.o
if [ $? -ne 0 ]; then
  echo "Errors while loading basis_compare.o"
  exit
fi
rm basis_compare.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/basis_compare
#
echo "Program installed as ~/bin/$ARCH/basis_compare"
