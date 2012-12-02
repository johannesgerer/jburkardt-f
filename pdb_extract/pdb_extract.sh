#!/bin/bash
#
gfortran -c -g pdb_extract.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pdb_extract.f90"
  exit
fi
rm compiler.txt
#
gfortran pdb_extract.o
if [ $? -ne 0 ]; then
  echo "Errors loading pdb_extract.o"
  exit
fi
#
rm pdb_extract.o
mv a.out ~/bin/$ARCH/pdb_extract
#
echo "Executable installed as ~/bin/$ARCH/pdb_extract"
