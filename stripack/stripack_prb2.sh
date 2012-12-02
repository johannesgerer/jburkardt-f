#!/bin/bash
#
gfortran -c -g stripack_prb2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stripack_prb2.f90"
  exit
fi
rm compiler.txt
#
gfortran stripack_prb2.o -L$HOME/lib/$ARCH -lstripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stripack_prb2.o"
  exit
fi
rm stripack_prb2.o
#
mv a.out stripack_prb2
./stripack_prb2 > stripack_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stripack_prb2"
  exit
fi
rm stripack_prb2
#
echo "Program output written to stripack_prb2_output.txt"
#
if [ -f stripack_prb2_del.eps ]; then
  convert stripack_prb2_del.eps stripack_prb2_del.png
  rm stripack_prb2_del.eps
fi
#
if [ -f stripack_prb2_vor.eps ]; then
  convert stripack_prb2_vor.eps stripack_prb2_vor.png
  rm stripack_prb2_vor.eps
fi
#
echo "Program graphics written to stripack_prb2_del.png and stripack_prb2_vor.png"
