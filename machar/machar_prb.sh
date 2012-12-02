#!/bin/bash
#
gfortran -c -g machar_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machar_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran machar_prb.o -L$HOME/lib/$ARCH -lmachar
if [ $? -ne 0 ]; then
  echo "Errors linking and loading machar_prb.o"
  exit
fi
rm machar_prb.o
#
mv a.out machar_prb
./machar_prb > machar_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running machar_prb"
  exit
fi
rm machar_prb
#
echo "Test program output written to machar_prb_output.txt."
