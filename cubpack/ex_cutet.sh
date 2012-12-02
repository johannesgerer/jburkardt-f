#!/bin/bash
#
ar x $HOME/lib/$ARCH/libcubpack.a precision_model.mod
ar x $HOME/lib/$ARCH/libcubpack.a cui.mod
#
gfortran -c -g ex_cutet.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ex_cutet.f90"
  exit
fi
rm compiler.txt
#
rm *.mod
#
gfortran ex_cutet.o -L$HOME/lib/$ARCH -lcubpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ex_cutet.o"
  exit
fi
rm ex_cutet.o
#
mv a.out ex_cutet
./ex_cutet > ex_cutet.out
if [ $? -ne 0 ]; then
  echo "Errors running ex_cutet"
  exit
fi
rm ex_cutet
#
echo "The ex_cutet test problem has been executed."
