#!/bin/bash
#
gfortran -c -g colored_noise_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling colored_noise_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran colored_noise_prb.o -L$HOME/lib/$ARCH -lcolored_noise
if [ $? -ne 0 ]; then
  echo "Errors linking and loading colored_noise_prb.o"
  exit
fi
rm colored_noise_prb.o
#
mv a.out colored_noise_prb
./colored_noise_prb > colored_noise_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running colored_noise_prb"
  exit
fi
rm colored_noise_prb
#
echo "Test program output written to colored_noise_prb_output.txt."
