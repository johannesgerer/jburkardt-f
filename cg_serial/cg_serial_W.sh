#!/bin/bash
#
cp npbparams_W.h npbparams.h
#
gfortran -c -g cg_serial.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cg_serial.f90"
  exit
fi
rm compiler.txt
#
gcc -c wtime.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime.c"
  exit
fi
rm compiler.txt
#
gfortran cg_serial.o wtime.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cg_serial.o + wtime.o"
  exit
fi
rm cg_serial.o
rm wtime.o
rm npbparams.h
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/cg_serial_W
#
echo "Executable installed as ~/bin/$ARCH/cg_serial_W"
