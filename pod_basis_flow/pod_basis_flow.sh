#!/bin/bash
#
gfortran -c -g pod_basis_flow.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pod_basis_flow.f90"
  exit
fi
rm compiler.txt
#
#gfortran pod_basis_flow.o -lcxml
gfortran pod_basis_flow.o -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pod_basis_flow.o"
  exit
fi
rm pod_basis_flow.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/pod_basis_flow
#
echo "Executable installed as ~/bin/$ARCH/pod_basis_flow"
