#!/bin/bash
#
gfortran -c -g snsimp.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling snsimp.f90"
  exit
fi
rm compiler.txt
#
gfortran snsimp.o -L$HOME/lib/$ARCH -larpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading snsimp.o"
  exit
fi
rm snsimp.o
#
mv a.out snsimp
./snsimp > snsimp_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running snsimp"
  exit
fi
rm snsimp
#
echo "Test program output written to snsimp_output.txt."
