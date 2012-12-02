#!/bin/bash
#
gfortran -c -g cvt_basis_flow.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_basis_flow.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_basis_flow.o -llapack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_basis_flow.o"
  exit
fi
rm cvt_basis_flow.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/cvt_basis_flow
#
echo "Program installed as ~/bin/$ARCH/cvt_basis_flow"
