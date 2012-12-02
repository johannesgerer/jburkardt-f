#!/bin/bash
#
gfortran -c -g codepack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling codepack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran codepack_prb.o -L$HOME/lib/$ARCH -lcodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading codepack_prb.o"
  exit
fi
rm codepack_prb.o
#
mv a.out codepack_prb
./codepack_prb > codepack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running codepack_prb"
  exit
fi
rm codepack_prb
#
echo "Test program output written to codepack_prb_output.txt."
