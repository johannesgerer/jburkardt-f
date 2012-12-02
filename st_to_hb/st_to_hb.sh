#!/bin/bash
#
gfortran -c -g st_to_hb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling st_to_hb.f90"
  exit
fi
rm compiler.txt
#
gfortran st_to_hb.o -L$HOME/lib/$ARCH -lhb_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading st_to_hb.o"
  exit
fi
rm st_to_hb.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/st_to_hb
#
echo "Executable installed as ~/bin/$ARCH/st_to_hb"
