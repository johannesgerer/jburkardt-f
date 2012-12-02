#!/bin/bash
#
gfortran -c -g diaphony.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling diaphony.f90"
  exit
fi
rm compiler.txt
#
gfortran diaphony.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading diaphony.o"
  exit
fi
#
rm diaphony.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/diaphony
#
echo "Executable installed as ~/bin/$ARCH/diaphony"
