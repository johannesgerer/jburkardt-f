#!/bin/bash
#
gfortran -c -g lattice_rule_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lattice_rule_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lattice_rule_prb.o -L$HOME/lib/$ARCH -llattice_rule
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lattice_rule_prb.o"
  exit
fi
rm lattice_rule_prb.o
#
mv a.out lattice_rule_prb
./lattice_rule_prb > lattice_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lattice_rule_prb"
  exit
fi
rm lattice_rule_prb
#
echo "Test program output written to lattice_rule_prb_output.txt."
