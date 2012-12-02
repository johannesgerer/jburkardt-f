#!/bin/bash
#
g95 -c -g g95_intrinsics.f90 >& compiler.txt
if [ $?-ne 0 ]; then
  echo "Errors compiling g95_intrinsics.f90"
  exit
fi
rm compiler.txt
#
g95 g95_intrinsics.o
if [ $?-ne 0 ]; then
  echo "Errors linking and loading g95_intrinsics.o"
  exit
fi
rm g95_intrinsics.o
#
mv a.out g95_intrinsics
./g95_intrinsics > g95_intrinsics_output.txt
if [ $?-ne 0 ]; then
  echo "Errors running g95_intrinsics"
  exit
fi
rm g95_intrinsics
#
echo "Program output written to g95_intrinsics_output.txt"
