#!/bin/bash
#
gfortran -c -g quadrature_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrature_test.f90"
  exit
fi
rm compiler.txt
#
gfortran quadrature_test.o -L$HOME/lib/$ARCH -ltest_nint
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrature_test.o + libtest_nint.a"
  exit
fi
rm quadrature_test.o
#
mv a.out $HOME/bin/$ARCH/quadrature_test
#
echo "Executable installed as ~/bin/$ARCH/quadrature_test"
