#!/bin/bash
#
gfortran -c pblas_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pblas_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran pblas_prb.o -L/$HOME/lib/$ARCH -lpblas -lblacs_f -lmpi
if [ $? -ne 0 ]; then
  echo "Errors compiling pblas_prb.o"
  exit
fi
rm *.o
#
mv a.out pblas_prb
pblas_prb > pblas_prb_output.txt
#
if [ $? -ne 0 ]; then
  echo "Errors running pblas_prb"
  exit
fi
rm pblas_prb
#
echo "Program output written to pblas_prb_output.txt"

