#!/bin/bash
#
gfortran -c -g ss_gd_align_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ss_gd_align_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ss_gd_align_prb.o -L$HOME/lib/$ARCH -lss_gd_align
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ss_gd_align_prb.o"
  exit
fi
rm ss_gd_align_prb.o
#
mv a.out ss_gd_align_prb
./ss_gd_align_prb > ss_gd_align_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ss_gd_align_prb"
  exit
fi
rm ss_gd_align_prb
#
echo "Test program output written to ss_gd_align_prb_output.txt."
