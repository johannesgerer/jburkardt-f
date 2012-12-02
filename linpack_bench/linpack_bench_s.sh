#!/bin/bash
#
gfortran -c -g linpack_bench_s.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_bench_s.f90"
  exit
fi
rm compiler.txt
#
gfortran linpack_bench_s.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_bench_s.o"
  exit
fi
rm linpack_bench_s.o
#
mv a.out ~/bin/$ARCH/linpack_bench_s
#
echo "Program installed as ~/bin/$ARCH/linpack_bench_s"
