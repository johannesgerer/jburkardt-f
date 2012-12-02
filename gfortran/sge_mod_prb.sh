#!/bin/bash
#
gfortran -c -g sge_mod.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sge_mod.f90"
  exit
fi
rm compiler.txt
#
gfortran -c -g sge_mod_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sge_mod_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sge_mod_prb.o sge_mod.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sge_mod_prb.o + sge_mod.o"
  exit
fi
rm sge_mod_prb.o
rm sge_mod.o
#
mv a.out sge_mod_prb
./sge_mod_prb > sge_mod_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sge_mod_prb"
  exit
fi
rm sge_mod_prb
#
echo "Program output written to sge_mod_prb_output.txt"
