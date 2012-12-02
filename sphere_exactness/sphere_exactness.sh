#!/bin/bash
#
gfortran -c -g sphere_exactness.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_exactness.f90"
  exit
fi
rm compiler.txt
#
gfortran sphere_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_exactness.o"
  exit
fi
rm sphere_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/sphere_exactness
#
echo "Program installed as ~/bin/$ARCH/sphere_exactness"
