#!/bin/bash
#
gfortran -c read_align.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling read_align.f90"
  exit
fi
rm compiler.txt
#
gfortran read_align.o
if [ $? -ne 0 ]; then
  echo "Errors loading read_align.o"
  exit
fi
#
rm read_align.o
mv a.out ~/bin/$ARCH/read_align
#
echo "Executable installed as ~/bin/$ARCH/read_align"
