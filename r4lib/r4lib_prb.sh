#!/bin/bash
#
gfortran -c -g r4lib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r4lib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran r4lib_prb.o -L$HOME/lib/$ARCH -lr4lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading r4lib_prb.o"
  exit
fi
rm r4lib_prb.o
#
mv a.out r4lib_prb
./r4lib_prb > r4lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running r4lib_prb"
  exit
fi
rm r4lib_prb
#
echo "Test program output written to r4lib_prb_output.txt."
