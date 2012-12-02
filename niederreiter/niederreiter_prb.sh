#!/bin/bash
#
gfortran -c -g niederreiter_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran niederreiter_prb.o -L$HOME/lib/$ARCH -lniederreiter
if [ $? -ne 0 ]; then
  echo "Errors linking and loading niederreiter_prb.o"
  exit
fi
rm niederreiter_prb.o
#
mv a.out niederreiter_prb
./niederreiter_prb > niederreiter_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running niederreiter_prb"
  exit
fi
rm niederreiter_prb
#
echo "Test program output written to niederreiter_prb_output.txt."
