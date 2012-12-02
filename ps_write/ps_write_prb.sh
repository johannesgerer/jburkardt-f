#!/bin/bash
#
gfortran -c -g ps_write_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ps_write_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ps_write_prb.o -L$HOME/lib/$ARCH -lps_write
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ps_write_prb.o"
  exit
fi
rm ps_write_prb.o
#
mv a.out ps_write_prb
./ps_write_prb > ps_write_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ps_write_prb"
  exit
fi
rm ps_write_prb
#
echo "Test program output written to ps_write_prb_output.txt."
