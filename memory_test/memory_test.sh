#!/bin/bash
#
gfortran -c -g memory_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling memory_test.f90"
  exit
fi
rm compiler.txt
#
gfortran memory_test.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading memory_test.o"
  exit
fi
rm memory_test.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/memory_test
#
echo "Program installed as ~/bin/$ARCH/memory_test"
