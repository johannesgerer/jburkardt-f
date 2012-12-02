#!/bin/bash
#
gfortran -c -g wandzura_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wandzura_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran wandzura_prb.o -L$HOME/lib/$ARCH -lwandzura
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wandzura_prb.o"
  exit
fi
rm wandzura_prb.o
#
mv a.out wandzura_prb
./wandzura_prb > wandzura_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wandzura_prb"
  exit
fi
rm wandzura_prb
#
echo "Test program output written to wandzura_prb_output.txt."
