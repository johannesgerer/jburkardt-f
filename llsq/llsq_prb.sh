#!/bin/bash
#
gfortran -c -g llsq_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling llsq_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran llsq_prb.o -L$HOME/lib/$ARCH -lllsq
if [ $? -ne 0 ]; then
  echo "Errors linking and loading llsq_prb.o"
  exit
fi
rm llsq_prb.o
#
mv a.out llsq_prb
./llsq_prb > llsq_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running llsq_prb"
  exit
fi
rm llsq_prb
#
echo "Test program output written to llsq_prb_output.txt."
