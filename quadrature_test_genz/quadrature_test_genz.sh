#!/bin/bash
#
gfortran -c -g quadrature_test_genz.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrature_test_genz.f90"
  exit
fi
rm compiler.txt
#
gfortran quadrature_test_genz.o -L$HOME/lib/$ARCH -ltest_nint
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrature_test_genz.o + libtest_nint.a"
  exit
fi
rm quadrature_test_genz.o
#
mv a.out $HOME/bin/$ARCH/quadrature_test_genz
#
echo "Executable installed as $HOME/bin/$ARCH/quadrature_test_genz"
