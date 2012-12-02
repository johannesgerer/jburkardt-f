#!/bin/bash
#
gfortran -c -g lake.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lake.f90"
  exit
fi
rm compiler.txt
#
gfortran ~/lib/$ARCH/fem2d_poisson_sparse.o lake.o -L$HOME/lib/$ARCH -lmgmres
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_sparse.o + lake.o"
  exit
fi
rm lake.o
#
chmod ugo+x a.out
mv a.out fem2d_poisson_sparse
fem2d_poisson_sparse lake > lake_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem2d_poisson_sparse."
  exit
fi
rm fem2d_poisson_sparse
#
if [ -e lake_elements.eps ]; then
  convert lake_elements.eps lake_elements.png
  rm lake_elements.eps
fi
#
if [ -e lake_nodes.eps ]; then
  convert lake_nodes.eps lake_nodes.png
  rm lake_nodes.eps
fi
#
echo "Program output written to lake_output.txt"
