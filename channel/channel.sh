#!/bin/bash
#
gfortran -c -g channel.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling channel.f90"
  exit
fi
rm compiler.txt
#
gfortran channel.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading channel.o"
  exit
fi
rm channel.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/channel
#
echo "Program installed as ~/bin/$ARCH/channel"
