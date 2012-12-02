#!/bin/bash
#
cp npbparams_S.h npbparams.h
#
gfortran -c -g ep_serial.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ep_serial.f90"
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
gfortran ep_serial.o wtime.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ep_serial.o + wtime.o"
  exit
fi
rm ep_serial.o
rm wtime.o
rm npbparams.h
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/ep_serial_S
#
echo "Executable installed as ~/bin/$ARCH/ep_serial_S"
