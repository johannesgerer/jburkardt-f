#!/bin/bash
#
gfortran -c -g sandia_sgmgg_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_sgmgg_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sandia_sgmgg_prb.o -L$HOME/lib/$ARCH -lsandia_sgmgg
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_sgmgg_prb.o"
  exit
fi
rm sandia_sgmgg_prb.o
#
mv a.out sandia_sgmgg_prb
./sandia_sgmgg_prb > sandia_sgmgg_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_sgmgg_prb"
  exit
fi
rm sandia_sgmgg_prb
#
echo "Test program output written to sandia_sgmgg_prb_output.txt."
