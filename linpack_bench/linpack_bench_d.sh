#!/bin/bash
#
gfortran -c -g linpack_bench_d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_bench_d.f90"
  exit
fi
rm compiler.txt
#
gfortran linpack_bench_d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_bench_d.o"
  exit
fi
rm linpack_bench_d.o
#
mv a.out ~/bin/$ARCH/linpack_bench_d
#
echo "Program installed as ~/bin/$ARCH/linpack_bench_d"
