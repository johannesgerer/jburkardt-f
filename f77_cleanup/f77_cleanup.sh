#!/bin/bash
#
gfortran -c -g f77_cleanup.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling f77_cleanup.f90"
  exit
fi
rm compiler.txt
#
gfortran f77_cleanup.o
if [ $? -ne 0 ]; then
  echo "Errors loading f77_cleanup.o"
  exit
fi
rm f77_cleanup.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/f77_cleanup
#
echo "Program installed as ~/bin/$ARCH/f77_cleanup."
