#!/bin/bash
#
gfortran -c anyplt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling anyplt_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran anyplt_prb.o -L$HOME/lib/$ARCH -lanynul
if [ $? -ne 0 ]; then
  echo "Errors linking and loading anyplt_prb.o"
  exit
fi
rm anyplt_prb.o
#
mv a.out anynul_prb
./anynul_prb > anynul_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running anynul_prb"
fi
rm anynul_prb
#
echo "Program output written to anynul_prb_output.txt"
