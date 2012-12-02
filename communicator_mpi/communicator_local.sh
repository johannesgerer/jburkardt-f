#!/bin/bash
#
mpif90 -c -g communicator_mpi.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling communicator_mpi.f90"
  exit
fi
rm compiler.txt
#
mpif90 communicator_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading communicator_mpi.o."
  exit
fi
#
rm communicator_mpi.o
#
mv a.out communicator
mpirun -np 4 ./communicator > communicator_local_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running communicator"
  exit
fi
rm communicator
#
echo "Program output written to communicator_local_output.txt"
