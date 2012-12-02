#!/bin/bash
#
gfortran -c lapack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lapack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lapack_prb.o -lscs
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lapack_prb.o"
  exit
fi
rm lapack_prb.o
#
mv a.out lapack_prb
./lapack_prb > lapack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lapack_prb"
  exit
fi
rm lapack_prb
#
echo "Program output written to lapack_prb_output.txt"
