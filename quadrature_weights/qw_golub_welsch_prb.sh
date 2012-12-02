#!/bin/bash
#
gfortran -c -g qw_golub_welsch_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qw_golub_welsch_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran qw_golub_welsch_prb.o -L$HOME/lib/$ARCH -lqw_golub_welsch
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qw_golub_welsch_prb.o"
  exit
fi
rm qw_golub_welsch_prb.o
#
mv a.out qw_golub_welsch_prb
./qw_golub_welsch_prb > qw_golub_welsch_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qw_golub_welsch_prb"
  exit
fi
rm qw_golub_welsch_prb
#
echo "Program output written to qw_golub_welsch_prb_output.txt"
