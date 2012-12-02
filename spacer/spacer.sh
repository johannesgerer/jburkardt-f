#!/bin/bash
#
gfortran -c spacer.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spacer.f90"
  exit
fi
rm compiler.txt
#
gfortran spacer.o
if [ $? -ne 0 ]; then
  echo "Errors linking spacer.o"
  exit
fi
#
rm spacer.o
mv a.out ~/bin/$ARCH/spacer
#
echo "Executable installed as ~/bin/$ARCH/spacer"
