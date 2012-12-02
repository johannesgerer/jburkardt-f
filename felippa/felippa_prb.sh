#!/bin/bash
#
gfortran -c -g felippa_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling felippa_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran felippa_prb.o -L$HOME/lib/$ARCH -lfelippa
if [ $? -ne 0 ]; then
  echo "Errors linking and loading felippa_prb.o"
  exit
fi
rm felippa_prb.o
#
mv a.out felippa_prb
./felippa_prb > felippa_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running felippa_prb"
  exit
fi
rm felippa_prb
#
echo "Test program output written to felippa_prb_output.txt."
