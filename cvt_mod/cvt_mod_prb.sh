#!/bin/bash
#
gfortran -c -g cvt_mod_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_mod_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_mod_prb.o -L$HOME/lib/$ARCH -lcvt_mod
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_mod_prb.o"
  exit
fi
rm cvt_mod_prb.o
#
mv a.out cvt_mod_prb
./cvt_mod_prb > cvt_mod_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_mod_prb"
  exit
fi
rm cvt_mod_prb
#
echo "Test program output written to cvt_mod_prb_output.txt."
