#!/bin/bash
#
mpif90 hello_mpi.f90
#
if [ $? -ne 0 ]; then
  echo "Errors compiling hello_mpi.f90"
  exit
fi
#
#  Rename the executable.
#
mv a.out hello
#
#  Ask MPI to use 8 processes to run your program.
#
mpirun -np 8 ./hello > hello_local_output.txt
#
#  Clean up.
#
rm hello
#
echo "Program output written to hello_local_output.txt"

