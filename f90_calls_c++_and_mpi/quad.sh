#!/bin/bash
#
source /usr/local/profile.d/openmpi-intel.sh
source /usr/local/profile.d/iccvars.sh
source /usr/local/profile.d/ifortvars.sh
#
mpif90 -c quad_main.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_main.f90"
  exit
fi
rm compiler.txt
#
mpiCC -c quad_sub.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_sub.C"
  exit
fi
rm compiler.txt
#
mpif90 quad_main.o quad_sub.o -lstdc++ -lmpi_cxx
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quad_main.o + quad_sub.o"
  exit
fi
#
mv a.out quad
#
./quad > quad_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quad."
  exit
fi
#
rm quad
#
echo "Program output written to quad_output.txt."
