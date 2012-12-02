#!/bin/bash
#
gfortran -c -g chebyshev1_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev1_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran chebyshev1_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev1_rule.o"
  exit
fi
rm chebyshev1_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/chebyshev1_rule
#
echo "Executable installed as ~/bin/$ARCH/chebyshev1_rule"
