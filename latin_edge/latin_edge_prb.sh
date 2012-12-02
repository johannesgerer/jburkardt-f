#!/bin/bash
#
gfortran -c -g latin_edge_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_edge_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran latin_edge_prb.o -L$HOME/lib/$ARCH -llatin_edge
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_edge_prb.o"
  exit
fi
rm latin_edge_prb.o
#
mv a.out latin_edge_prb
./latin_edge_prb > latin_edge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latin_edge_prb"
  exit
fi
rm latin_edge_prb
#
echo "Test program output written to latin_edge_prb_output.txt."
