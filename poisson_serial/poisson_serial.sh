#!/bin/bash
#
gfortran -c poisson_serial.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_serial.f90"
  exit
fi
rm compiler.txt
#
gfortran poisson_serial.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading poisson_serial.o"
  exit
fi
rm poisson_serial.o
#
mv a.out poisson_serial
poisson_serial > poisson_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poisson_serial"
  exit
fi
rm poisson_serial
#
echo "Program output written to poisson_serial_output.txt"
