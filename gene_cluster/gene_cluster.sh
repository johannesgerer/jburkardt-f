#!/bin/bash
#
gfortran -c -g gene_cluster.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling gene_cluster.f90"
  exit
fi
rm compiler.txt
#
gfortran gene_cluster.o
if [ $? -ne 0 ]; then
  echo "Errors while loading gene_cluster.o"
  exit
fi
rm gene_cluster.o
#
mv a.out ~/bin/$ARCH/gene_cluster
#
echo "Program installed as ~/bin/$ARCH/gene_cluster"
