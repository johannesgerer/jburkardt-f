#!/bin/bash
#
gfortran -c f77_to_f90.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling f77_to_f90.f90"
  exit
fi
rm compiler.txt
#
gfortran f77_to_f90.o
if [ $? -ne 0 ]; then
  echo "Errors while loading f77_to_f90.o"
  exit
fi
rm f77_to_f90.o
rm *.mod
#
mv a.out ~/bin/$ARCH/f77_to_f90
#
echo "Program installed as ~/bin/$ARCH/f77_to_f90"
