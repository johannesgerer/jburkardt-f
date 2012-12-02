#!/bin/bash
#
ar x $HOME/lib/$ARCH/libcubpack.a precision_model.mod
ar x $HOME/lib/$ARCH/libcubpack.a cui.mod
#
gfortran -c -g ex_triex.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ex_triex.f90"
  exit
fi
rm compiler.txt
#
rm *.mod
#
gfortran ex_triex.o -L$HOME/lib/$ARCH -lcubpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ex_triex.o"
  exit
fi
rm ex_triex.o
#
mv a.out ex_triex
./ex_triex > ex_triex.out
if [ $? -ne 0 ]; then
  echo "Errors running ex_triex"
  exit
fi
rm ex_triex
#
echo "The ex_triex test problem has been executed."
