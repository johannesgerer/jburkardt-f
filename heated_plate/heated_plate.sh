#!/bin/bash
#
gfortran -c -g heated_plate.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling heated_plate.f90"
  exit
fi
rm compiler.txt
#
gfortran heated_plate.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading heated_plate.o"
  exit
fi
rm heated_plate.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/heated_plate
#
echo "The program has been installed as ~/bin/$ARCH/heated_plate."
