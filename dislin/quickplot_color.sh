#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran quickplot_color.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quickplot_color.f90"
  exit
fi
rm compiler.txt
#
gfortran quickplot_color.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quickplot_color.o"
  exit
fi
rm quickplot_color.o
#
mv a.out quickplot_color
./quickplot_color > quickplot_color_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quickplot_color"
  exit
fi
rm quickplot_color
#
echo "Program output written to quickplot_color_output.txt"
