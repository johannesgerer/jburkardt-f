#!/bin/bash
#
gfortran ../hello_openmp/hello_openmp.f90 -L$HOME/lib/$ARCH -lopenmp_stubs
mv a.out hello
#
#  Run the program.
#
./hello > hello_output.txt
rm hello
#
echo "Program output written to hello_output.txt"
