#!/bin/bash
#
gfortran -c -g test_triangulation_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_triangulation_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_triangulation_prb.o -L$HOME/lib/$ARCH -ltest_triangulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_triangulation_prb.o"
  exit
fi
rm test_triangulation_prb.o
#
mv a.out test_triangulation_prb
./test_triangulation_prb > test_triangulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_triangulation_prb"
  exit
fi
rm test_triangulation_prb
#
echo "Test program output written to test_triangulation_prb_output.txt."
