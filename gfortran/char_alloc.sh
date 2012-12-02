#!/bin/bash
#
gfortran -c -g char_alloc.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling char_alloc.f90"
  exit
fi
rm compiler.txt
#
gfortran char_alloc.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading char_alloc.o"
  exit
fi
rm char_alloc.o
#
mv a.out char_alloc
./char_alloc > char_alloc_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running char_alloc"
  exit
fi
rm char_alloc
#
echo "Program output written to char_alloc_output.txt"
