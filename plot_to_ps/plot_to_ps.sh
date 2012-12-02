#!/bin/bash
#
gfortran -c -g plot_to_ps.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling plot_to_ps.f90"
  exit
fi
rm compiler.txt
#
gfortran plot_to_ps.o
if [ $? -ne 0 ]; then
  echo "Errors linking plot_to_ps.o"
  exit
fi
#
rm plot_to_ps.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/plot_to_ps
#
echo "Program installed as ~/bin/$ARCH/plot_to_ps"
