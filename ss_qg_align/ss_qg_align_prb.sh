#!/bin/bash
#
gfortran -c -g ss_qg_align_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ss_qg_align_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ss_qg_align_prb.o -L$HOME/lib/$ARCH -lss_qg_align
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ss_qg_align_prb.o"
  exit
fi
rm ss_qg_align_prb.o
#
mv a.out ss_qg_align_prb
./ss_qg_align_prb > ss_qg_align_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ss_qg_align_prb"
  exit
fi
rm ss_qg_align_prb
#
echo "Test program output written to ss_qg_align_prb_output.txt."
