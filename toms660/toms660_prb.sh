#!/bin/bash
#
gfortran -c -g toms660_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms660_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms660_prb.o -L$HOME/lib/$ARCH -ltoms660
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms660_prb.o"
  exit
fi
rm toms660_prb.o
#
mv a.out toms660_prb
./toms660_prb > toms660_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms660_prb"
  exit
fi
rm toms660_prb
#
echo "Test program output written to toms660_prb_output.txt."
