#!/bin/bash
#
gfortran -c file_merge.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while compiling file_merge.f90"
  exit
fi
rm compiler.txt
#
gfortran file_merge.o
if [ $? -ne 0 ]; then
  echo "Errors occurred while linking file_merge.o"
  exit
fi
rm file_merge.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/file_merge
#
echo "Program installed as ~/bin/$ARCH/file_merge"
