#!/bin/bash
#
gfortran -c -g ziggurat_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ziggurat_prb.o -L$HOME/lib/$ARCH -lziggurat
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ziggurat_prb.o"
  exit
fi
rm ziggurat_prb.o
#
mv a.out ziggurat_prb
./ziggurat_prb > ziggurat_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ziggurat_prb"
  exit
fi
rm ziggurat_prb
#
echo "Test program output written to ziggurat_prb_output.txt."
