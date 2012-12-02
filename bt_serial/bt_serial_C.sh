#!/bin/bash
#
cp npbparams_C.h npbparams.h
#
gfortran -c -g bt_serial.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bt_serial.f90"
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
gfortran bt_serial.o wtime.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bt_serial.o + wtime.o"
  exit
fi
rm bt_serial.o
rm wtime.o
rm npbparams.h
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/bt_serial_C
#
echo "Executable installed as ~/bin/$ARCH/bt_serial_C"
