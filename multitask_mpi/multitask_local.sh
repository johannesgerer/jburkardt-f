#!/bin/bash
#
mpif90 multitask_mpi.f90
#
if [ $? -ne 0 ]; then
  echo "Errors compiling multitask_mpi.f90"
  exit
fi
#
#  Rename the executable.
#
mv a.out multitask
#
#  Ask MPI to use 8 processes to run your program.
#
mpirun -np 8 ./multitask > multitask_local_output.txt
#
#  Clean up.
#
rm multitask
#
echo "Program output written to multitask_local_output.txt"

