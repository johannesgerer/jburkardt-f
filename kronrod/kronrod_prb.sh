#!/bin/bash
#
gfortran -c -g kronrod_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kronrod_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran kronrod_prb.o -L$HOME/lib/$ARCH -lkronrod
if [ $? -ne 0 ]; then
  echo "Errors linking and loading kronrod_prb.o"
  exit
fi
rm kronrod_prb.o
#
mv a.out kronrod_prb
./kronrod_prb > kronrod_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running kronrod_prb"
  exit
fi
rm kronrod_prb
#
echo "Test program output written to kronrod_prb_output.txt."
