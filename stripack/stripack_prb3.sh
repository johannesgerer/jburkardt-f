#!/bin/bash
#
gfortran -c -g stripack_prb3.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stripack_prb3.f90"
  exit
fi
rm compiler.txt
#
gfortran stripack_prb3.o -L$HOME/lib/$ARCH -lstripack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stripack_prb3.o"
  exit
fi
rm stripack_prb3.o
#
mv a.out stripack_prb3
./stripack_prb3 > stripack_prb3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stripack_prb3"
  exit
fi
rm stripack_prb3
#
echo "Program output written to stripack_prb3_output.txt"
#
if [ -f stripack_prb3_del.eps ]; then
  convert stripack_prb3_del.eps stripack_prb3_del.png
  rm stripack_prb3_del.eps
fi
#
if [ -f stripack_prb3_vor.eps ]; then
  convert stripack_prb3_vor.eps stripack_prb3_vor.png
  rm stripack_prb3_vor.eps
fi
#
echo "Program graphics written to stripack_prb3_del.png and stripack_prb3_vor.png"
