#!/bin/bash
#
gfortran -c -g ps_gg_align_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ps_gg_align_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ps_gg_align_prb.o -L$HOME/lib/$ARCH -lps_gg_align
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ps_gg_align_prb.o"
  exit
fi
rm ps_gg_align_prb.o
#
mv a.out ps_gg_align_prb
./ps_gg_align_prb > ps_gg_align_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ps_gg_align_prb"
  exit
fi
rm ps_gg_align_prb
#
echo "Test program output written to ps_gg_align_prb_output.txt."
