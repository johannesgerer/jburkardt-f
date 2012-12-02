#!/bin/bash
#
#  Compile with GFORTRAN.
#
gfortran -fopenmp dot_product.f90
#
#  Compile with IFORT.
#
#ifort -openmp -parallel -fpp dot_product.f90
#
mv a.out dot_product
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./dot_product > dot_product_local_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./dot_product >> dot_product_local_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./dot_product >> dot_product_local_output.txt
#
#  Discard the executable.
#
rm dot_product
#
echo "Program output written to dot_product_local_output.txt"