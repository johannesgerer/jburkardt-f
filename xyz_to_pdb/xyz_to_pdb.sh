#!/bin/bash
#
gfortran -c -g xyz_to_pdb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling xyz_to_pdb.f90"
  exit
fi
rm compiler.txt
#
gfortran xyz_to_pdb.o
if [ $? -ne 0 ]; then
  echo "Error loading xyz_to_pdb.o"
  exit
fi
rm xyz_to_pdb.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/xyz_to_pdb
#
echo "Executable installed as ~/bin/$ARCH/xyz_to_pdb"
