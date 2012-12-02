#!/bin/bash
#
gfortran -c -g fem3d_pack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem3d_pack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran fem3d_pack_prb.o -L$HOME/lib/$ARCH -lfem3d_pack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem3d_pack_prb.o"
  exit
fi
rm fem3d_pack_prb.o
#
mv a.out fem3d_pack_prb
./fem3d_pack_prb > fem3d_pack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem3d_pack_prb"
  exit
fi
rm fem3d_pack_prb
#
echo "Test program output written to fem3d_pack_prb_output.txt."
