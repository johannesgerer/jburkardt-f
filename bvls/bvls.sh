#!/bin/bash
#
gfortran -c -g bvls.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bvls_prb.f90"
  exit
fi
rm compiler.txt
#
mv bvls.o ~/lib/$ARCH
#
echo "Object code installed as ~/lib/$ARCH/bvls.o"
