#!/bin/bash
#
gfortran -c -g distance_to_position.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling distance_to_position.f90"
  exit
fi
rm compiler.txt
#
gfortran distance_to_position.o -L$HOME/lib/$ARCH -lnms
if [ $? -ne 0 ]; then
  echo "Errors linking and loading distance_to_position.o"
  exit
fi
rm distance_to_position.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/distance_to_position
#
echo "The executable has been installed as ~/bin/$ARCH/distance_to_position"
