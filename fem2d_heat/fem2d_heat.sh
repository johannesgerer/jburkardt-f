#!/bin/bash
#
gfortran -c -g fem2d_heat.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_heat.f90"
  exit
fi
rm compiler.txt
#
mv fem2d_heat.o ~/lib/$ARCH
#
echo "Partial program object code installed as ~/lib/$ARCH/fem2d_heat.o"
