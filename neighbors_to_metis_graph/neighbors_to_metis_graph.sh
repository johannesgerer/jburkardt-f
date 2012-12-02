#!/bin/bash
#
gfortran -c -g neighbors_to_metis_graph.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling neighbors_to_metis_graph.f90"
  exit
fi
rm compiler.txt
#
gfortran neighbors_to_metis_graph.o
if [ $? -ne 0 ]; then
  echo "Error loading neighbors_to_metis_graph.o
  exit
fi
rm neighbors_to_metis_graph.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/neighbors_to_metis_graph
#
echo "Executable installed as ~/bin/$ARCH/neighbors_to_metis_graph"
