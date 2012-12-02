#!/bin/bash
#
gfortran -c -g mandelbrot.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot.f90"
  exit
fi
rm compiler.txt
#
gfortran mandelbrot.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mandelbrot.o"
  exit
fi
rm mandelbrot.o
#
mv a.out ~/bin/$ARCH/mandelbrot
#
echo "Executable installed as ~/bin/$ARCH/mandelbrot"
