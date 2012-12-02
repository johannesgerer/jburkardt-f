#!/bin/bash
#
gfortran -c xyz_plot.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling xyz_plot.f90"
  exit
fi
#
gfortran xyz_plot.o ~/lib/$ARCH/libanynul.a
if [ $? -ne 0 ]; then
  echo "Errors loading xyz_plot.o + libanynul.a"
  exit
fi
rm xyz_plot.o
#
mv a.out ~/bin/$ARCH/xyz_plot_nul
#
echo "Executable installed as ~/bin/$ARCH/xyz_plot_nul"
