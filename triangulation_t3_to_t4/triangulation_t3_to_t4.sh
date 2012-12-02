#!/bin/bash
#
gfortran -c -g triangulation_t3_to_t4.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_t3_to_t4.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_t3_to_t4.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_t3_to_t4.o"
  exit
fi
#
rm triangulation_t3_to_t4.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_t3_to_t4
#
echo "Executable installed as ~/bin/$ARCH/triangulation_t3_to_t4"
