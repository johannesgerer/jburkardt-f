#!/bin/bash
#
gfortran -c -g tripack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tripack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran tripack_prb.o -L$HOME/lib/$ARCH -ltripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tripack_prb.o"
  exit
fi
rm tripack_prb.o
#
mv a.out tripack_prb
./tripack_prb > tripack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tripack_prb"
  exit
fi
rm tripack_prb
#
echo "Test program output written to tripack_prb_output.txt."
