#!/bin/bash
#
g95 -c -g real_kind16.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling real_kind16.f90"
  exit
fi
rm compiler.txt
#
g95 real_kind16.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading real_kind16.o"
  exit
fi
rm real_kind16.o
#
mv a.out real_kind16
./real_kind16 > real_kind16_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running real_kind16"
  exit
fi
rm real_kind16
#
echo "The real_kind16 test problem has been executed."
