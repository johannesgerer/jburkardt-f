#!/bin/bash
#
g95 -c -g g95_quadmath.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling g95_quadmath.f90"
  exit
fi
rm compiler.txt
#
g95 g95_quadmath.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading g95_quadmath.o"
  exit
fi
rm g95_quadmath.o
#
mv a.out g95_quadmath
./g95_quadmath > g95_quadmath_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running g95_quadmath"
# exit
fi
rm g95_quadmath
#
echo "Program output written to g95_quadmath_output.txt."
