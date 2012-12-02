#!/bin/bash
#
gfortran -c -g files_multiple.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling files_multiple.f90"
  exit
fi
rm compiler.txt
#
gfortran files_multiple.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading files_multiple.o"
  exit
fi
rm files_multiple.o
#
mv a.out files_multiple
./files_multiple > files_multiple_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running files_multiple"
  exit
fi
rm files_multiple
#
echo "Program output written to files_multiple_output.txt."
