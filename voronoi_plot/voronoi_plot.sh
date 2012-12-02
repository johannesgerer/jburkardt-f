#!/bin/bash
#
gfortran -c -g voronoi_plot.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling voronoi_plot.f90"
  exit
fi
rm compiler.txt
#
gfortran voronoi_plot.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading voronoi_plot.o"
  exit
fi
rm voronoi_plot.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/voronoi_plot
#
echo "Executable installed as ~/bin/$ARCH/voronoi_plot"
