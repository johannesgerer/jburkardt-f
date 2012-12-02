#!/bin/bash
#
gfortran -c -g integral_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling integral_test.f90"
  exit
fi
rm compiler.txt
#
gfortran integral_test.o -L$HOME/lib/$ARCH -ltest_nint
if [ $? -ne 0 ]; then
  echo "Errors linking and loading integral_test.o"
  exit
fi
rm integral_test.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/integral_test
#
echo "Executable installed as ~/bin/$ARCH/integral_test"
