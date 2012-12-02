#!/bin/bash
#
gfortran -c -g cvt_mod_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_mod_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_mod_dataset.o -L$HOME/lib/$ARCH -lcvt_mod
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_mod_dataset.o + libcvt_mod.a"
  exit
fi
rm cvt_mod_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/cvt_mod_dataset
#
echo "Program installed as ~/bin/$ARCH/cvt_mod_dataset"
