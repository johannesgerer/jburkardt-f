#!/bin/bash
#
gfortran -c -g getwgt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling getwgt_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran getwgt_prb.o -L$HOME/lib/$ARCH -lgetwgt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading getwgt_prb.o"
  exit
fi
rm getwgt_prb.o
#
mv a.out getwgt_prb
./getwgt_prb > getwgt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running getwgt_prb"
  exit
fi
rm getwgt_prb
#
echo "Test program output written to getwgt_prb_output.txt."
