#!/bin/bash
#
gfortran -c scalapack_prb1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling scalapack_prb1.f90"
  exit
fi
rm compiler.txt
#
gfortran scalapack_prb1.o -L$HOME/lib/$ARCH -lscalapack -lblacs_f
if [ $? -ne 0 ]; then
  echo "Errors linking and loading scalapack_prb1.o"
  exit
fi
rm scalapack_prb1.o
#
mv a.out scalapack_prb1
./scalapack_prb1 > scalapack_prb1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running scalapack_prb1"
  exit
fi
rm scalapack_prb1
#
echo "Program output written to scalapack_prb1_output.txt"
