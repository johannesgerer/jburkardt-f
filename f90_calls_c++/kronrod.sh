#!/bin/bash
#
gfortran -c -g kronrod_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kronrod_prb.f90"
  exit
fi
rm compiler.txt
#
g++ -c -g kronrod.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kronrod_prb.cpp"
  exit
fi
rm compiler.txt
#
gfortran kronrod_prb.o kronrod.o -lm -lstdc++
if [ $? -ne 0 ]; then
  echo "Errors linking and loading kronrod_prb.o + kronrod.o"
  exit
fi
rm kronrod_prb.o
rm kronrod.o
#
mv a.out kronrod
./kronrod > kronrod_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running kronrod"
  exit
fi
rm kronrod
#
echo "Test program output written to kronrod_output.txt."
