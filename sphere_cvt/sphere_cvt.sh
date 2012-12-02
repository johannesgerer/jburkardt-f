#!/bin/bash
#
gfortran -c -g sphere_cvt.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_cvt.f90"
  exit
fi
rm compiler.txt
#
gfortran sphere_cvt.o -L$HOME/lib/$ARCH/ -lstripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_cvt.o + stripack.a"
  exit
fi
rm sphere_cvt.o
#
mv a.out ~/bin/$ARCH/sphere_cvt
#
echo "Executable installed as ~/bin/$ARCH/sphere_cvt"
