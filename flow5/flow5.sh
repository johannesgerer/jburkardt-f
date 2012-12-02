#!/bin/bash
#
gfortran -c flow5.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling flow5.f90"
  exit
fi
rm compiler.txt
#
gfortran flow5.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading flow5.o"
  exit
fi
rm flow5.o
#
mv a.out ~/bin/$ARCH/flow5
#
echo "Program installed as ~/bin/$ARCH/flow5"
