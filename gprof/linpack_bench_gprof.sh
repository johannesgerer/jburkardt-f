#!/bin/bash
#
gfortran -c -pg linpack_bench.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_bench.f90"
  exit
fi
rm compiler.txt
#
gfortran -pg linpack_bench.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_bench.o"
  exit
fi
rm linpack_bench.o
#
mv a.out linpack_bench
#
./linpack_bench > linpack_bench_output.txt
#
gprof linpack_bench >> linpack_bench_gprof_output.txt
#
#  Clean up.
#  GPROF creates a temporary file GMON.OUT that we don't need.
#
rm linpack_bench
rm gmon.out
#
echo "Program output written to linpack_bench_output.txt"
echo "GPROF report written to linpack_bench_gprof_output.txt"
