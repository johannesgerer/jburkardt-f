#!/bin/bash
#
gfortran -c -g sandia_sparse_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_sparse_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sandia_sparse_prb.o -L$HOME/lib/$ARCH -lsandia_sparse
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_sparse_prb.o"
  exit
fi
rm sandia_sparse_prb.o
#
mv a.out sandia_sparse_prb
./sandia_sparse_prb > sandia_sparse_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_sparse_prb"
  exit
fi
rm sandia_sparse_prb
#
echo "Test program output written to sandia_sparse_prb_output.txt."
