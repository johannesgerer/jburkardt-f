#!/bin/bash
#
gfortran -c -g ps_qg_align_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ps_qg_align_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ps_qg_align_prb.o -L$HOME/lib/$ARCH -lps_qg_align
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ps_qg_align_prb.o"
  exit
fi
rm ps_qg_align_prb.o
#
mv a.out ps_qg_align_prb
./ps_qg_align_prb > ps_qg_align_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ps_qg_align_prb"
  exit
fi
rm ps_qg_align_prb
#
echo "Test program output written to ps_qg_align_prb_output.txt."
