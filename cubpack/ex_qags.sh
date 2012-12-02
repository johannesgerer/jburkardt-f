#!/bin/bash
#
ar x $HOME/lib/$ARCH/libcubpack.a precision_model.mod
ar x $HOME/lib/$ARCH/libcubpack.a cui.mod
#
gfortran -c -g ex_qags.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ex_qags.f90"
  exit
fi
rm compiler.txt
#
rm *.mod
#
gfortran ex_qags.o -L$HOME/lib/$ARCH -lcubpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ex_qags.o"
  exit
fi
rm ex_qags.o
#
mv a.out ex_qags
./ex_qags > ex_qags.out
if [ $? -ne 0 ]; then
  echo "Errors running ex_qags"
  exit
fi
rm ex_qags
#
echo "The ex_qags test problem has been executed."
