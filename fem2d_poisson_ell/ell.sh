#!/bin/bash
#
gfortran -c -g ell.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ell.f90"
  exit
fi
rm compiler.txt
#
gfortran ~/lib/$ARCH/fem2d_poisson.o ell.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson.o + ell.o"
  exit
fi
rm ell.o
#
chmod ugo+x a.out
mv a.out fem2d_poisson
fem2d_poisson ell > ell_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem2d_poisson."
  exit
fi
rm fem2d_poisson
#
if [ -e ell_elements.eps ]; then
  convert ell_elements.eps ell_elements.png
  rm ell_elements.eps
fi
#
if [ -e ell_nodes.eps ]; then
  convert ell_nodes.eps ell_nodes.png
  rm ell_nodes.eps
fi
#
echo "Program output written to ell_output.txt"
