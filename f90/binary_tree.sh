#!/bin/bash
#
gfortran -c -g binary_tree.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling binary_tree.f90"
  exit
fi
rm compiler.txt
#
gfortran binary_tree.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading binary_tree.o"
  exit
fi
rm binary_tree.o
rm *.mod
#
mv a.out binary_tree
./binary_tree > binary_tree_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running binary_tree"
  exit
fi
rm binary_tree
#
echo "Program output written to binary_tree_output.txt"
