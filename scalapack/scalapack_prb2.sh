#!/bin/bash
#
gfortran -c scalapack_prb2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling scalapack_prb2.f90"
  exit
fi
rm compiler.txt
#
gfortran scalapack_prb2.o -L$HOME/lib/$ARCH -lscalapack -lblacs_f
if [ $? -ne 0 ]; then
  echo "Errors linking and loading scalapack_prb2.o"
  exit
fi
rm scalapack_prb2.o
#
mv a.out scalapack_prb2
./scalapack_prb2 > scalapack_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running scalapack_prb2"
  exit
fi
rm scalapack_prb2
#
echo "Program output written to scalapack_prb2_output.txt"
