#!/bin/bash
#
gfortran -c -g square.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling square.f90"
  exit
fi
rm compiler.txt
#
#  Link the precompiled main program FEM2D_HEAT with the user routines.
#
gfortran ~/lib/$ARCH/fem2d_heat.o square.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_heat.o + square.o"
  exit
fi
rm square.o
#
#  Run the program with the user mesh files.
#
chmod ugo+x a.out
mv a.out fem2d_heat_square
fem2d_heat_square square_nodes.txt square_elements.txt > square_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem2d_heat_square."
  exit
fi
rm fem2d_heat_square
#
#  Convert the EPS files to PNG format.
#
if [ -e elements.eps ]; then
  convert elements.eps square_elements.png
  rm elements.eps
fi
#
if [ -e nodes.eps ]; then
  convert nodes.eps square_nodes.png
  rm nodes.eps
fi
#
echo "Program output written to square_output.txt"
