#!/bin/bash
#
gfortran -c -g include_files.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling include_files.f90"
  exit
fi
rm compiler.txt
#
gfortran include_files.o
if [ $? -ne 0 ]; then
  echo "Errors loading include_files.o"
  exit
fi
rm include_files.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/include_files
#
echo "Executable installed as ~/bin/$ARCH/include_files"
