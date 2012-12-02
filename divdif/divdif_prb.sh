#!/bin/bash
#
gfortran -c -g divdif_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling divdif_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran divdif_prb.o -L$HOME/lib/$ARCH -ldivdif
if [ $? -ne 0 ]; then
  echo "Errors linking and loading divdif_prb.o"
  exit
fi
rm divdif_prb.o
#
mv a.out divdif_prb
./divdif_prb > divdif_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running divdif_prb"
  exit
fi
rm divdif_prb
#
echo "Program output written to divdif_prb_output.txt"
