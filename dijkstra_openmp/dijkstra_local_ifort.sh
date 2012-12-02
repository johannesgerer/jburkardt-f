#!/bin/bash
#
#  Compile with IFORT.
#
ifort -openmp -parallel -fpp dijkstra_openmp.f90
#
mv a.out dijkstra
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./dijkstra > dijkstra_local_ifort_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./dijkstra >> dijkstra_local_ifort_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./dijkstra >> dijkstra_local_ifort_output.txt
#
#  Discard the executable.
#
rm dijkstra
#
echo "Program output written to dijkstra_local_ifort_output.txt"
