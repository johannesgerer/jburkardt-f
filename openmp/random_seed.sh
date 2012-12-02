#!/bin/bash
#
gfortran -fseed random_seed.f90
mv a.out random_seed
#
./random_seed > random_seed_output.txt
rm random_seed
#
echo "Program output written to random_seed_output.txt."
