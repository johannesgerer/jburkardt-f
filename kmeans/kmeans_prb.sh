#!/bin/bash
#
gfortran -c -g kmeans_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kmeans_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran kmeans_prb.o -L$HOME/lib/$ARCH -lkmeans
if [ $? -ne 0 ]; then
  echo "Errors linking and loading kmeans_prb.o"
  exit
fi
rm kmeans_prb.o
#
mv a.out kmeans_prb
./kmeans_prb > kmeans_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running kmeans_prb"
  exit
fi
rm kmeans_prb
#
echo "Test program output written to kmeans_prb_output.txt."
