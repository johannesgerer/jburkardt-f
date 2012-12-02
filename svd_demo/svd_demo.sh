#!/bin/bash
#
gfortran -c -g svd_demo.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_demo.f90"
  exit
fi
rm compiler.txt
#
gfortran svd_demo.o -L$HOME/lib/$ARCH -llapack -llinpack_d -lblas1_d
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_demo.o"
  exit
fi
rm svd_demo.o
#
mv a.out svd_demo
mv svd_demo $HOME/bin/$ARCH
#
echo "Executable installed as ~/bin/$ARCH/svd_demo"
