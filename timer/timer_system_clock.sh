#!/bin/bash
#
gfortran -c -g timer_system_clock.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timer_system_clock.f90"
  exit
fi
rm compiler.txt
#
gfortran timer_system_clock.o
if [ $? -ne 0 ]; then
  echo "Errors loading timer_system_clock.o"
  exit
fi
rm timer_system_clock.o
#
mv a.out timer_system_clock
./timer_system_clock > timer_system_clock_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running timer_system_clock"
  exit
fi
rm timer_system_clock
#
echo "Program output written to timer_system_clock_output.txt."
