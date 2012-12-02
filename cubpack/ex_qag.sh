#!/bin/bash
#
ar x $HOME/lib/$ARCH/libcubpack.a precision_model.mod
ar x $HOME/lib/$ARCH/libcubpack.a cui.mod
#
gfortran -c -g ex_qag.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ex_qag.f90"
  exit
fi
rm compiler.txt
#
rm *.mod
#
gfortran ex_qag.o -L$HOME/lib/$ARCH -lcubpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ex_qag.o"
  exit
fi
rm ex_qag.o
#
mv a.out ex_qag
./ex_qag > ex_qag.out
if [ $? -ne 0 ]; then
  echo "Errors running ex_qag"
  exit
fi
rm ex_qag
#
echo "The ex_qag test problem has been executed."
