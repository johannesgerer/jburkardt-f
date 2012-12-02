#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran scatter_plot.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling scatter_plot.f90"
  exit
fi
rm compiler.txt
#
gfortran scatter_plot.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading scatter_plot.o"
  exit
fi
rm scatter_plot.o
#
mv a.out scatter_plot
./scatter_plot > scatter_plot_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running scatter_plot"
  exit
fi
rm scatter_plot
#
echo "Program output written to scatter_plot_output.txt"
