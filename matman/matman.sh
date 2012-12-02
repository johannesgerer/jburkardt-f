#!/bin/bash
#
gfortran -g -c matman.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matman.f90"
  exit
fi
rm compiler.txt
#
gfortran matman.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matman.o"
  exit
fi
rm matman.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/matman
#
echo "Program installed as ~/bin/$ARCH/matman"
