#!/bin/bash
#
gfortran -c -g owens_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling owens_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran owens_prb.o -L$HOME/lib/$ARCH -lowens
if [ $? -ne 0 ]; then
  echo "Errors linking and loading owens_prb.o"
  exit
fi
rm owens_prb.o
#
mv a.out owens_prb
./owens_prb > owens_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running owens_prb"
  exit
fi
rm owens_prb
#
echo "Test program output written to owens_prb_output.txt."
