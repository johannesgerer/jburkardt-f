#!/bin/bash
#
gfortran -c -g tetrahedron_monte_carlo_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_monte_carlo_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran tetrahedron_monte_carlo_prb.o -L$HOME/lib/$ARCH -ltetrahedron_monte_carlo
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_monte_carlo_prb.o"
  exit
fi
rm tetrahedron_monte_carlo_prb.o
#
mv a.out tetrahedron_monte_carlo_prb
./tetrahedron_monte_carlo_prb > tetrahedron_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tetrahedron_monte_carlo_prb"
  exit
fi
rm tetrahedron_monte_carlo_prb
#
echo "Test program output written to tetrahedron_monte_carlo_prb_output.txt."
