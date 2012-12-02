#!/bin/bash
#
gfortran -c -g vector_plot.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vector_plot.f90"
  exit
fi
rm compiler.txt
#
gfortran vector_plot.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vector_plot.o"
  exit
fi
#
rm vector_plot.o
mv a.out ~/bin/$ARCH/vector_plot
#
echo "Executable installed as ~/bin/$ARCH/vector_plot"
