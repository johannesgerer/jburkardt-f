#!/bin/bash
#
gfortran -c -g vec_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vec_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran vec_io_prb.o -L$HOME/lib/$ARCH -lvec_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vec_io_prb.o"
  exit
fi
rm vec_io_prb.o
#
mv a.out vec_io_prb
./vec_io_prb > vec_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running vec_io_prb"
  exit
fi
rm vec_io_prb
#
echo "Test program output written to vec_io_prb_output.txt."
