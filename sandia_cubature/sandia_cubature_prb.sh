#!/bin/bash
#
gfortran -c -g sandia_cubature_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_cubature_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sandia_cubature_prb.o -L$HOME/lib/$ARCH -lsandia_cubature
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_cubature_prb.o"
  exit
fi
rm sandia_cubature_prb.o
#
mv a.out sandia_cubature_prb
./sandia_cubature_prb > sandia_cubature_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_cubature_prb"
  exit
fi
rm sandia_cubature_prb
#
echo "Test program output written to sandia_cubature_prb_output.txt."
