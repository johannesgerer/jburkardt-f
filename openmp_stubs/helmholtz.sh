#!/bin/bash
#
gfortran ../openmp/helmholtz.f90 -L$HOME/lib/$ARCH -lopenmp_stubs
mv a.out helmholtz
#
#  Run the program.
#
./helmholtz > helmholtz_output.txt
rm helmholtz
#
echo "Program output written to helmholtz_output.txt"
