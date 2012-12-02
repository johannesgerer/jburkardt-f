#!/bin/bash
#
gfortran -c -g bdmlib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bdmlib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran bdmlib_prb.o -L$HOME/lib/$ARCH -lbdmlib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bdmlib_prb.o"
  exit
fi
rm bdmlib_prb.o
#
mv a.out bdmlib_prb
./bdmlib_prb > bdmlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bdmlib_prb"
  exit
fi
rm bdmlib_prb
#
echo "Test program output written to bdmlib_prb_output.txt."
