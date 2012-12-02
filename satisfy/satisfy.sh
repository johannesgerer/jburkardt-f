#!/bin/bash
#
gfortran -c -g satisfy.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling satisfy.f90."
  exit
fi
rm compiler.txt
#
gfortran satisfy.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading satisfy.o."
  exit
fi
#
rm satisfy.o
#
mv a.out satisfy
./satisfy > satisfy_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running satisfy."
  exit
fi
rm satisfy
#
echo "Program output written to satisfy_output.txt"
