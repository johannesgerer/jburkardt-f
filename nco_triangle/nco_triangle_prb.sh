#!/bin/bash
#
gfortran -c -g nco_triangle_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nco_triangle_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran nco_triangle_prb.o -L$HOME/lib/$ARCH -lnco_triangle
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nco_triangle_prb.o"
  exit
fi
rm nco_triangle_prb.o
#
mv a.out nco_triangle_prb
./nco_triangle_prb > nco_triangle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nco_triangle_prb"
  exit
fi
rm nco_triangle_prb
#
echo "Test program output written to nco_triangle_prb_output.txt."
