#!/bin/bash
#
gfortran -c -g plot3d_to_avs.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling plot3d_to_avs.f90"
  exit
fi
rm compiler.txt
#
gfortran plot3d_to_avs.o -L$HOME/lib/$ARCH -lplot3d_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading plot3d_to_avs.o"
  exit
fi
rm plot3d_to_avs.o
#
mv a.out ~/bin/$ARCH/plot3d_to_avs
echo "Library installed as ~/bin/$ARCH/plot3d_to_avs"
