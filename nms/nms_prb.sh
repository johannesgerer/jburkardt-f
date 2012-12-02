#!/bin/bash
#
gfortran -c -g nms_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nms_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran nms_prb.o -L$HOME/lib/$ARCH -lnms
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nms_prb.o"
  exit
fi
rm nms_prb.o
#
mv a.out nms_prb
./nms_prb > nms_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nms_prb"
  exit
fi
rm nms_prb
#
echo "Test program output written to nms_prb_output.txt."
