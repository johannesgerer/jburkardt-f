#!/bin/bash
#
gfortran -c -g ncc_triangle_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ncc_triangle_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ncc_triangle_prb.o -L$HOME/lib/$ARCH -lncc_triangle
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ncc_triangle_prb.o"
  exit
fi
rm ncc_triangle_prb.o
#
mv a.out ncc_triangle_prb
./ncc_triangle_prb > ncc_triangle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ncc_triangle_prb"
  exit
fi
rm ncc_triangle_prb
#
echo "Test program output written to ncc_triangle_prb_output.txt."
