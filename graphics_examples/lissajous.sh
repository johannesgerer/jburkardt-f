#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran lissajous.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lissajous.f90"
  exit
fi
rm compiler.txt
#
gfortran lissajous.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lissajous.o"
  exit
fi
rm lissajous.o
#
mv a.out lissajous
./lissajous > lissajous_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lissajous"
  exit
fi
rm lissajous
#
echo "Program output written to lissajous_output.txt"
