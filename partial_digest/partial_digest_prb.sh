#!/bin/bash
#
gfortran -c -g partial_digest_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling partial_digest_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran partial_digest_prb.o -L$HOME/lib/$ARCH -lpartial_digest
if [ $? -ne 0 ]; then
  echo "Errors linking and loading partial_digest_prb.o"
  exit
fi
rm partial_digest_prb.o
#
mv a.out partial_digest_prb
./partial_digest_prb > partial_digest_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running partial_digest_prb"
  exit
fi
rm partial_digest_prb
#
echo "Test program output written to partial_digest_prb_output.txt."
