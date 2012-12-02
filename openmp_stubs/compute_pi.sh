#!/bin/bash
#
gfortran ../openmp/compute_pi.f90 -L$HOME/lib/$ARCH -lopenmp_stubs
mv a.out compute_pi
#
#  Run the program.
#
./compute_pi > compute_pi_output.txt
rm compute_pi
#
echo "Program output written to compute_pi_output.txt"
