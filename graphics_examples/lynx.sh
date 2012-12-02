#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran lynx.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lynx.f90"
  exit
fi
rm compiler.txt
#
gfortran lynx.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lynx.o"
  exit
fi
rm lynx.o
#
mv a.out lynx
./lynx > lynx_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lynx"
  exit
fi
rm lynx
#
echo "Program output written to lynx_output.txt"
