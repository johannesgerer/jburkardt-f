#!/bin/bash
#
gfortran -c -g file_name_sequence_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling file_name_sequence_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran file_name_sequence_prb.o -L$HOME/lib/$ARCH -lfile_name_sequence
if [ $? -ne 0 ]; then
  echo "Errors linking and loading file_name_sequence_prb.o"
  exit
fi
rm file_name_sequence_prb.o
#
mv a.out file_name_sequence_prb
./file_name_sequence_prb > file_name_sequence_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running file_name_sequence_prb"
  exit
fi
rm file_name_sequence_prb
#
echo "Test program output written to file_name_sequence_prb_output.txt."
