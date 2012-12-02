#!/bin/bash
#
gfortran -c -g linked_list.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linked_list.f90"
  exit
fi
rm compiler.txt
#
gfortran linked_list.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linked_list.o"
  exit
fi
rm linked_list.o
rm *.mod
#
mv a.out linked_list
./linked_list > linked_list_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linked_list"
  exit
fi
rm linked_list
#
echo "Program output written to linked_list_output.txt"
