#!/bin/bash
#
gfortran -c -g pdb_read_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pdb_read_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran pdb_read_prb.o -L$HOME/lib/$ARCH -lpdb_read
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pdb_read_prb.o"
  exit
fi
rm pdb_read_prb.o
#
mv a.out pdb_read_prb
./pdb_read_prb < pdb_read_prb_input.txt > pdb_read_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pdb_read_prb"
  exit
fi
rm pdb_read_prb
#
echo "Test program output written to pdb_read_prb_output.txt."
