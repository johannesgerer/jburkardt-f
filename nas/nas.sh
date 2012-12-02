#!/bin/bash
#
gfortran -c -g nas.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nas.f90"
  exit
fi
rm compiler.txt
#
gfortran nas.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nas.o"
  exit
fi
rm nas.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/nas
#
echo "Executable installed as ~/bin/$ARCH/nas"
