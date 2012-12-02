#!/bin/bash
#
gfortran -c -g calpak_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling calpak_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran calpak_prb.o -L$HOME/lib/$ARCH -lcalpak
if [ $? -ne 0 ]; then
  echo "Errors linking and loading calpak_prb.o"
  exit
fi
rm calpak_prb.o
#
mv a.out calpak_prb
./calpak_prb > calpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running calpak_prb"
  exit
fi
rm calpak_prb
#
echo "Test program output written to calpak_prb_output.txt."
