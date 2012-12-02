#!/bin/bash
#
#  Compile the program with IFORT.
#
ifort -openmp -parallel -fpp ziggurat_openmp.f90
#
mv a.out ziggurat
#
#  Run with 1, 2, 4 and 8 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./ziggurat > ziggurat_local_ifort_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./ziggurat >> ziggurat_local_ifort_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./ziggurat >> ziggurat_local_ifort_output.txt
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./ziggurat >> ziggurat_local_ifort_output.txt
#
#  Discard the executable file.
#
rm ziggurat
#
echo "Program output written to ziggurat_local_ifort_output.txt"
