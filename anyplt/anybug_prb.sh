#!/bin/bash
#
gfortran -c -g anyplt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran anyplt_prb.o -L$HOME/lib/$ARCH -lanybug
if [ $? -ne 0 ]; then
  echo "Errors linking and loading anyplt_prb.o"
  exit
fi
rm anyplt_prb.o
#
mv a.out anybug_prb
./anybug_prb > anybug_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running anybug_prb"
fi
rm anybug_prb
mv anyplt.bug anybug_prb.bug
#
echo "Program output written to anybug_prb_output.txt"
