#!/bin/bash
#
gfortran -c -g toms659_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms659_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms659_prb.o -L$HOME/lib/$ARCH -ltoms659
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms659_prb.o"
  exit
fi
rm toms659_prb.o
#
mv a.out toms659_prb
./toms659_prb < toms659_prb_input.txt > toms659_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms659_prb"
  exit
fi
rm toms659_prb
#
echo "Test program output written to toms659_prb_output.txt."
