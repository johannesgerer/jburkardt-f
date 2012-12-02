#!/bin/bash
#
gfortran -c -g latin_edge_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_edge_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran latin_edge_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_edge_dataset.o + liblatin_edge.a"
  exit
fi
rm latin_edge_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/latin_edge_dataset
#
echo "Executable installed as ~/bin/$ARCH/latin_edge_dataset"
