#!/bin/bash
#
gfortran -c -g sweep2_voronoi_eps.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sweep2_voronoi_eps.f90"
  exit
fi
rm compiler.txt
#
gfortran sweep2_voronoi_eps.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sweep2_voronoi_eps.o"
  exit
fi
rm sweep2_voronoi_eps.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/sweep2_voronoi_eps
#
echo "Executable installed as ~/bin/$ARCH/sweep2_voronoi_eps"
