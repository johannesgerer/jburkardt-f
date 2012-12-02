#!/bin/bash
#
gfortran -c md3glue.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md3glue.f90"
  exit
fi
rm compiler.txt
#
gfortran md3glue.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md3glue.o"
  exit
fi
rm md3glue.o
#
mv a.out ~/bin/$ARCH/md3glue
#
echo "Executable installed as ~/bin/$ARCH/md3glue"
