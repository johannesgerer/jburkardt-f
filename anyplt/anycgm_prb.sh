#!/bin/bash
#
gfortran -c -g anyplt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran anyplt_prb.o -L$HOME/lib/$ARCH -lanycgm -ldrawcgm -lXt -lX11 -lXmu
if [ $? -ne 0 ]; then
  echo "Errors linking and loading anyplt_prb.o"
  exit
fi
rm anyplt_prb.o
#
mv a.out anycgm_prb
./anycgm_prb > anycgm_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running anycgm_prb"
fi
rm anycgm_prb
#
echo "Program output written to anycgm_prb_output.txt"
