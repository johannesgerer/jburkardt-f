#!/bin/bash
#
gfortran -c -g fem2d_project_function.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_project_function.f90"
  exit
fi
rm compiler.txt
#
gfortran fem2d_project_function.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_project_function.o"
  exit
fi
rm fem2d_project_function.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/fem2d_project_function
#
echo "Program installed as ~/bin/$ARCH/fem2d_project_function"
