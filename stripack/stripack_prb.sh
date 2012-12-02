#!/bin/bash
#
gfortran -c -g stripack_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stripack_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran stripack_prb.o -L$HOME/lib/$ARCH -lstripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stripack_prb.o"
  exit
fi
rm stripack_prb.o
#
mv a.out stripack_prb
./stripack_prb > stripack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stripack_prb"
  exit
fi
rm stripack_prb
#
echo "Program output written to stripack_prb_output.txt"
#
if [ -f stripack_prb_del.eps ]; then
  convert stripack_prb_del.eps stripack_prb_del.png
  rm stripack_prb_del.eps
fi
#
if [ -f stripack_prb_vor.eps ]; then
  convert stripack_prb_vor.eps stripack_prb_vor.png
  rm stripack_prb_vor.eps
fi
#
echo "Program graphics written to stripack_prb_del.png and stripack_prb_vor.png"
