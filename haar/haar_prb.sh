#!/bin/bash
#
gfortran -c -g haar_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling haar_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran haar_prb.o -L$HOME/lib/$ARCH -lhaar
if [ $? -ne 0 ]; then
  echo "Errors linking and loading haar_prb.o"
  exit
fi
rm haar_prb.o
#
mv a.out haar_prb
./haar_prb > haar_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running haar_prb"
  exit
fi
rm haar_prb
#
echo "Test results written to haar_prb_output.txt."
