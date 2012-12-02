#!/bin/bash
#
ar x $HOME/lib/$ARCH/libcubpack.a precision_model.mod
ar x $HOME/lib/$ARCH/libcubpack.a cui.mod
#
gfortran -c -g details.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling details.f90"
  exit
fi
rm compiler.txt
#
rm *.mod
#
gfortran details.o -L$HOME/lib/$ARCH -lcubpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading details.o"
  exit
fi
rm details.o
#
mv a.out details
./details > details.out
if [ $? -ne 0 ]; then
  echo "Errors running details"
  exit
fi
rm details
#
echo "The details test problem has been executed."
