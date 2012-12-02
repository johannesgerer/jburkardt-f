#!/bin/bash
#
gfortran -c -g bar_plot_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bar_plot_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran bar_plot_prb.o -L$HOME/lib/$ARCH -lbar_plot
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bar_plot_prb.o"
  exit
fi
rm bar_plot_prb.o
#
mv a.out bar_plot_prb
./bar_plot_prb > bar_plot_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bar_plot_prb"
  exit
fi
rm bar_plot_prb
#
echo "Test program output written to bar_plot_prb_output.txt."
