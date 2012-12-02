#!/bin/bash
#
cp ../../datasets/tsp/att48.tsp .
cp ../../datasets/tsp/p01.tsp .
#
gfortran -c -g tsp_io_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tsp_io_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran tsp_io_prb.o -L$HOME/lib/$ARCH -ltsp_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tsp_io_prb.o"
  exit
fi
rm tsp_io_prb.o
#
mv a.out tsp_io_prb
./tsp_io_prb > tsp_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tsp_io_prb"
  exit
fi
#
mv p01_d.txt ../../datasets/tsp
mv att48_d.txt ../../datasets/tsp
#
rm tsp_io_prb
rm att48.tsp
rm p01.tsp
#
echo "Test program output written to tsp_io_prb_output.txt."
