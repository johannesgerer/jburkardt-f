#!/bin/bash
#
gfortran -c -g quad2d_serial.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad2d_serial.f90"
  exit
fi
rm compiler.txt
#
gfortran quad2d_serial.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quad2d_serial.o"
  exit
fi
rm quad2d_serial.o
#
mv a.out quad2d_serial
./quad2d_serial > quad2d_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quad2d_serial"
  exit
fi
rm quad2d_serial
#
echo "Program output written to quad2d_serial_output.txt"
