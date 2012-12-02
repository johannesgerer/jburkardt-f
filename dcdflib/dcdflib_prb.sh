#!/bin/bash
#
gfortran -c -g dcdflib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dcdflib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran dcdflib_prb.o -L$HOME/lib/$ARCH -ldcdflib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dcdflib_prb.o"
  exit
fi
rm dcdflib_prb.o
#
mv a.out dcdflib_prb
./dcdflib_prb > dcdflib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dcdflib_prb"
  exit
fi
rm dcdflib_prb
#
echo "Test program output written to dcdflib_prb_output.txt."
