#!/bin/bash
#
gfortran -c -g quadpack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadpack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran quadpack_prb.o -L$HOME/lib/$ARCH -lquadpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadpack_prb.o"
  exit
fi
rm quadpack_prb.o
#
mv a.out quadpack_prb
./quadpack_prb > quadpack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadpack_prb"
  exit
fi
rm quadpack_prb
#
echo "Test program output written to quadpack_prb_output.txt."
