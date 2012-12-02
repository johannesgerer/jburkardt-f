#!/bin/bash
#
mpif90 search_mpi.f90
#
if [ $? -ne 0 ]; then
  echo "Errors compiling search_mpi.f90"
  exit
fi
#
#  Rename the executable.
#
mv a.out search
#
#  Ask MPI to use 8 processes to run your program.
#
mpirun -np 8 ./search > search_local_output.txt
#
#  Clean up.
#
rm search
#
echo "Program output written to search_local_output.txt"

