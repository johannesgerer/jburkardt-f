#!/bin/bash
#
gfortran -c arby4.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling arby4.f90"
  exit
fi
rm compiler.txt
#
gfortran arby4.o -L$HOME/lib/$ARCH -ltoms611
if [ $? -ne 0 ]; then
  echo "Errors linking and loading arby4.o"
  exit
fi
rm arby4.o
#
mv a.out ~/bin/$ARCH/arby4
#
echo "Executable stored as ~/bin/$ARCH/arby4"
