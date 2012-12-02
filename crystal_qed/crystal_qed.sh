#!/bin/bash
#
gfortran -c crystal_qed.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling crystal_qed.f90"
  exit
fi
rm compiler.txt
#
gfortran crystal_qed.o -L$HOME/lib/$ARCH -ldqed
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading crystal_qed.o"
  exit
fi
rm crystal_qed.o
#
mv a.out ~/bin/$ARCH/crystal_qed
#
echo "Executable installed as ~/bin/$ARCH/crystal_qed"
