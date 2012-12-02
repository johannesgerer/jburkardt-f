#!/bin/bash
#
gfortran -c -g quadrature_test_2d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrature_test_2d.f90"
  exit
fi
rm compiler.txt
#
gfortran quadrature_test_2d.o -L$HOME/lib/$ARCH -ltest_int_2d
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrature_test_2d.o + libtest_int_2d.a"
  exit
fi
rm quadrature_test_2d.o
#
mv a.out $HOME/bin/$ARCH/quadrature_test_2d
#
echo "Executable installed as ~/bin/$ARCH/quadrature_test_2d"
