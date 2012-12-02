#!/bin/bash
#
gfortran -c -g steam_interact.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling steam_interact.f90"
  exit
fi
rm compiler.txt
#
gfortran steam_interact.o -L$HOME/lib/$ARCH -lsteam
if [ $? -ne 0 ]; then
  echo "Errors linking and loading steam_interact.o"
  exit
fi
rm steam_interact.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/steam_interact
#
echo "Executable installed as ~/bin/$ARCH/steam_interact"
