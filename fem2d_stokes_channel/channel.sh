#!/bin/bash
#
gfortran -c -g channel.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling channel.f90"
  exit
fi
rm compiler.txt
#
gfortran ~/lib/$ARCH/fem2d_stokes.o channel.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_stokes.o + channel.o"
  exit
fi
rm channel.o
#
chmod ugo+x a.out
mv a.out channel
./channel nodes6.txt triangles6.txt > channel_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running channel."
  exit
fi
rm channel
#
echo "Program output written to channel_output.txt"
