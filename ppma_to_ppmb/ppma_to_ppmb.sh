#!/bin/bash
#
gfortran -c -g ppma_to_ppmb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppma_to_ppmb.f90"
  exit
fi
rm compiler.txt
#
gfortran ppma_to_ppmb.o
if [ $? -ne 0 ]; then
  echo "Errors loading ppma_to_ppmb.o"
  exit
fi
#
rm ppma_to_ppmb.o
mv a.out ~/bin/$ARCH/ppma_to_ppmb
#
echo "Executable installed as ~/bin/$ARCH/ppma_to_ppmb"
