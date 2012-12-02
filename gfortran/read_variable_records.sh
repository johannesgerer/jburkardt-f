#!/bin/bash
#
gfortran -c -g read_variable_records.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling read_variable_records.f90"
  exit
fi
rm compiler.txt
#
gfortran read_variable_records.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading read_variable_records.o"
  exit
fi
rm read_variable_records.o
#
mv a.out read_variable_records
./read_variable_records > read_variable_records_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running read_variable_records"
  exit
fi
rm read_variable_records
#
echo "Program output written to read_variable_records_output.txt"
