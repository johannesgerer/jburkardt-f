#!/bin/bash
#
gfortran -c -g ihs_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ihs_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran ihs_dataset.o -L$HOME/lib/$ARCH -lihs
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ihs_dataset.o + libihs.a"
  exit
fi
rm ihs_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/ihs_dataset
#
echo "Executable installed as ~/bin/$ARCH/ihs_dataset"
