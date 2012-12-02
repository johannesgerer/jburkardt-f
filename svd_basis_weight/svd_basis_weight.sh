#!/bin/bash
#
gfortran -c -g svd_basis_weight.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_basis_weight.f90"
  exit
fi
rm compiler.txt
#
#gfortran svd_basis_weight.o -lcxml
gfortran svd_basis_weight.o -framework veclib
#gfortran svd_basis_weight.o -llapack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_basis_weight.o"
  exit
fi
rm svd_basis_weight.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/svd_basis_weight
#
echo "Executable installed as ~/bin/$ARCH/svd_basis_weight"
