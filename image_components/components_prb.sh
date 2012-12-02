#!/bin/bash
#
gfortran -c -g components_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling components_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran components_prb.o -L$HOME/lib/$ARCH -lcomponents
if [ $? -ne 0 ]; then
  echo "Errors linking and loading components_prb.o"
  exit
fi
rm components_prb.o
#
mv a.out components_prb
./components_prb > components_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running components_prb"
  exit
fi
rm components_prb
#
echo "Test program output written to components_prb_output.txt."
