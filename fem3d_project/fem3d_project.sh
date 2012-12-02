#!/bin/bash
#
gfortran -c fem3d_project.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem3d_project.f90"
  exit
fi
rm compiler.txt
#
gfortran fem3d_project.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem3d_project.o"
  exit
fi
#
rm fem3d_project.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/fem3d_project
#
echo "Program installed as ~/bin/$ARCH/fem3d_project"
