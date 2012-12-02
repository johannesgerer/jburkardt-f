#!/bin/bash
#
gfortran -c -g jbar.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jbar.f90"
  exit
fi
rm compiler.txt
#
gfortran jbar.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading jbar.o"
  exit
fi
rm jbar.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/jbar
#
echo "Executable installed as ~/bin/$ARCH/jbar"
