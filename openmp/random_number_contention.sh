#!/bin/bash
#
gfortran -fopenmp random_number_contention.f90
mv a.out random_number_contention
#
./random_number_contention > random_number_contention_output.txt
rm random_number_contention
#
echo "Program output written to random_number_contention_output.txt"
