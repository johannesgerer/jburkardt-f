#!/bin/bash
#
gfortran -c -g shepard_interp_1d_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shepard_interp_1d_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran shepard_interp_1d_prb.o -L$HOME/lib/$ARCH -lshepard_interp_1d -ltest_interp -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shepard_interp_1d_prb.o"
  exit
fi
rm shepard_interp_1d_prb.o
#
mv a.out shepard_interp_1d_prb
./shepard_interp_1d_prb > shepard_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running shepard_interp_1d_prb"
  exit
fi
rm shepard_interp_1d_prb
#
echo "Test program output written to shepard_interp_1d_prb_output.txt."
