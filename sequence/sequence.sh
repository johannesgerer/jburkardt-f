#!/bin/bash
#
gfortran -c sequence.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sequence.f90"
  exit
fi
rm compiler.txt
#
gfortran sequence.o
if [ $? -ne 0 ]; then
  echo "Errors loading sequence.o"
  exit
fi
rm sequence.o
#
mv a.out ~/bin/$ARCH/sequence
#
echo "Program installed as ~/bin/$ARCH/sequence"
