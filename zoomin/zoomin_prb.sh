#!/bin/bash
#
gfortran -c -g zoomin_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling zoomin_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran zoomin_prb.o -L$HOME/lib/$ARCH -lzoomin
if [ $? -ne 0 ]; then
  echo "Errors linking and loading zoomin_prb.o"
  exit
fi
rm zoomin_prb.o
#
mv a.out zoomin_prb
./zoomin_prb > zoomin_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running zoomin_prb"
  exit
fi
rm zoomin_prb
#
echo "Test program output written to zoomin_prb_output.txt."
