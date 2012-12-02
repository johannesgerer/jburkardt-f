#!/bin/bash
#
gfortran -c flow3.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while compiling flow3.f90"
  exit
fi
rm compiler.txt
#
gfortran flow3.o -L$HOME/lib/$ARCH -ltoms611
if [ $? -ne 0 ]; then
  echo "Errors occurred while linking and loading flow3.o"
  exit
fi
rm flow3.o
#
mv a.out ~/bin/$ARCH/flow3
#
echo "Program installed as ~/bin/$ARCH/flow3"
