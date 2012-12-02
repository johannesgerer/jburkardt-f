#!/bin/bash
#
XLF90 -c -g xlf_intrinsics_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xlf_intrinsics_prb.f90"
  exit
fi
rm compiler.txt
#
XLF90 xlf_intrinsics_prb.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xlf_intrinsics_prb.o"
  exit
fi
rm xlf_intrinsics_prb.o
#
mv a.out xlf_intrinsics_prb
./xlf_intrinsics_prb > xlf_intrinsics_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running xlf_intrinsics_prb"
  exit
fi
rm xlf_intrinsics_prb
#
echo "Program output written to xlf_intrinsics_prb_output.txt"
