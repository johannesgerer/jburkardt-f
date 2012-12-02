#!/bin/bash
#
gfortran -c -g fem2d_stokes.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_stokes.f90"
  exit
fi
rm compiler.txt
#
mv fem2d_stokes.o ~/lib/$ARCH
#
echo "Object code installed as ~/lib/$ARCH/fem2d_stokes.o"
