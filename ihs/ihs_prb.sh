#!/bin/bash
#
gfortran -c -g ihs_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ihs_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran ihs_prb.o -L$HOME/lib/$ARCH -lihs
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ihs_prb.o"
  exit
fi
rm ihs_prb.o
#
mv a.out ihs_prb
./ihs_prb > ihs_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ihs_prb"
  exit
fi
rm ihs_prb
#
echo "Test program output written to ihs_prb_output.txt."
