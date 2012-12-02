#!/bin/bash
#
gfortran -c -g plot_points.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling plot_points.f90"
  exit
fi
rm compiler.txt
#
gfortran -g plot_points.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading plot_points.o"
  exit
fi
#
rm plot_points.o
#
#  Because you are using a module, you're likely to have a funny
#  compiled module file lying around.
#
rm *.mod
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/plot_points
#
echo "Executable installed as ~/bin/$ARCH/plot_points"
