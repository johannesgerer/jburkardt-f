#!/bin/bash
#
gfortran -c -g timer_cpu_time.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timer_cpu_time.f90"
  exit
fi
rm compiler.txt
#
gfortran timer_cpu_time.o
if [ $? -ne 0 ]; then
  echo "Errors loading timer_cpu_time.o"
  exit
fi
rm timer_cpu_time.o
#
mv a.out timer_cpu_time
./timer_cpu_time > timer_cpu_time_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running timer_cpu_time"
  exit
fi
rm timer_cpu_time
#
echo "Program output written to timer_cpu_time_output.txt"
