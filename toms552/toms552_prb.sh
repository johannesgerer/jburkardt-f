#!/bin/bash
#
gfortran -c -g toms552_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms552_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms552_prb.o -L$HOME/lib/$ARCH -ltoms552
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms552_prb.o"
  exit
fi
rm toms552_prb.o
#
mv a.out toms552_prb
./toms552_prb < toms552_prb_input.txt > toms552_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms552_prb"
  exit
fi
rm toms552_prb
#
echo "Test program output written to toms552_prb_output.txt."
