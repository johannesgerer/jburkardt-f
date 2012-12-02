#!/bin/bash
#
gfortran -c -g tanh_sinh_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tanh_sinh_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran tanh_sinh_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tanh_sinh_rule.o"
  exit
fi
rm tanh_sinh_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/tanh_sinh_rule
#
echo "Executable installed as ~/bin/$ARCH/tanh_sinh_rule"
