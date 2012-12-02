#!/bin/bash
#
gfortran -c -g toms708_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms708_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms708_prb.o -L$HOME/lib/$ARCH -ltoms708
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms708_prb.o"
  exit
fi
rm toms708_prb.o
#
mv a.out toms708_prb
./toms708_prb > toms708_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms708_prb"
  exit
fi
rm toms708_prb
#
echo "Test program output written to toms708_prb_output.txt."
