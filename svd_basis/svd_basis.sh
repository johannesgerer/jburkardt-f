#!/bin/bash
#
gfortran -c -g svd_basis.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_basis.f90"
  exit
fi
rm compiler.txt
#
#gfortran svd_basis.o -L$HOME/lib/$ARCH -llapack
#gfortran svd_basis.o -L$HOME/lib/$ARCH -lcxml
gfortran svd_basis.o -L$HOME/lib/$ARCH -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_basis.o"
  exit
fi
rm svd_basis.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/svd_basis
#
echo "Executable installed as ~/bin/$ARCH/svd_basis"
