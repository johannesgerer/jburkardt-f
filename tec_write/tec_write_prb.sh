#!/bin/bash
#
gfortran -c -g tec_write_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tec_write_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran tec_write_prb.o -L$HOME/lib/$ARCH -ltec_write
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tec_write_prb.o"
  exit
fi
rm tec_write_prb.o
#
mv a.out tec_write_prb
./tec_write_prb > tec_write_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tec_write_prb"
  exit
fi
rm tec_write_prb
#
echo "Test program output written to tec_write_prb_output.txt."
