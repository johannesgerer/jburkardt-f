#!/bin/bash
#
gfortran -c -g crystal_plot.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  exit
fi
rm compiler.txt
#
gfortran crystal_plot.o ~/lib/$ARCH/libdrawcgm.a -L/usr/X11/lib -lXt -lX11 -lXmu
if [ $? -ne 0 ]; then
  exit
fi
#
rm crystal_plot.o
mv a.out ~/bin/$ARCH/crystal_plot
#
echo "Executable installed as ~/bin/$ARCH/crystal_plot."
