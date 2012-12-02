#!/bin/bash
#
gfortran -c -g sparsepak_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsepak_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsepak_prb.o -L$HOME/lib/$ARCH -lsparsepak
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsepak_prb.o"
  exit
fi
rm sparsepak_prb.o
#
mv a.out sparsepak_prb
./sparsepak_prb > sparsepak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsepak_prb"
  exit
fi
rm sparsepak_prb
#
echo "Test program output written to sparsepak_prb_output.txt."
