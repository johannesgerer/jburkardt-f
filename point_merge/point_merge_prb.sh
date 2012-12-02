#!/bin/bash
#
gfortran -c -g point_merge_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling point_merge_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran point_merge_prb.o -L$HOME/lib/$ARCH -lpoint_merge
if [ $? -ne 0 ]; then
  echo "Errors linking and loading point_merge_prb.o"
  exit
fi
rm point_merge_prb.o
#
mv a.out point_merge_prb
./point_merge_prb > point_merge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running point_merge_prb"
  exit
fi
rm point_merge_prb
#
echo "Test program output written to point_merge_prb_output.txt."
