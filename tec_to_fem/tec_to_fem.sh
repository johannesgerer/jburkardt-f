#!/bin/bash
#
gfortran -c -g tec_to_fem.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tec_to_fem.f90"
  exit
fi
rm compiler.txt
#
gfortran tec_to_fem.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tec_to_fem.o"
  exit
fi
#
rm tec_to_fem.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/tec_to_fem
#
echo "Executable installed as ~/bin/$ARCH/tec_to_fem"
