#!/bin/bash
#
gfortran -c -g stripack_interactive.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stripack_interactive.f90"
  exit
fi
rm compiler.txt
#
gfortran stripack_interactive.o -L$HOME/lib/$ARCH -lstripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stripack_interactive.o"
  exit
fi
rm stripack_interactive.o
#
chmod ugo+x a.out
mv a.out $HOME/bin/$ARCH/stripack_interactive
#
echo "Executable installed as ~/bin/$ARCH/stripack_interactive"
