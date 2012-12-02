#!/bin/bash
#
gfortran -c matvec.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matvec.f90"
  exit
fi
rm compiler.txt
#
gfortran matvec.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matvec.o"
  exit
fi
rm matvec.o
#
mv a.out matvec
mpirun -v -np 4 ./matvec > matvec_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running matvec"
  exit
fi
rm matvec
#
echo "Program output written to matvec_output.txt"
