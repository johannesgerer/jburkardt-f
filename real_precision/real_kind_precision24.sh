#!/bin/bash
#
g95 -c -g real_kind_precision24.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling real_kind_precision24.f90"
  exit
fi
rm compiler.txt
#
g95 real_kind_precision24.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading real_kind_precision24.o"
  exit
fi
rm real_kind_precision24.o
#
mv a.out real_kind_precision24
./real_kind_precision24 > real_kind_precision24_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running real_kind_precision24"
  exit
fi
rm real_kind_precision24
#
echo "The real_kind_precision24 test problem has been executed."
