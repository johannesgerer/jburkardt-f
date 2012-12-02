#!/bin/bash
#
gfortran -c module_mark.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling module_mark.f90"
  exit
fi
rm compiler.txt
#
gfortran module_mark.o
if [ $? -ne 0 ]; then
  echo "Errors loading module_mark.o"
  exit
fi
rm module_mark.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/module_mark
#
echo "Executable installed as ~/bin/$ARCH/module_mark"
