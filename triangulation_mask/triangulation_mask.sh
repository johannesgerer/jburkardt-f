#!/bin/bash
#
gfortran -c -g triangulation_mask.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_mask.f90"
  exit
fi
rm compiler.txt
#
mv triangulation_mask.o ~/lib/$ARCH
#
echo "Partial program installed as ~/lib/$ARCH/triangulation_mask.o"
