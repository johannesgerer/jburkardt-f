#!/bin/bash
#
gfortran -c -g lyness_rule_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lyness_rule_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lyness_rule_prb.o -L$HOME/lib/$ARCH -llyness_rule
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lyness_rule_prb.o"
  exit
fi
rm lyness_rule_prb.o
#
mv a.out lyness_rule_prb
./lyness_rule_prb > lyness_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lyness_rule_prb"
  exit
fi
rm lyness_rule_prb
#
echo "Test program output written to lyness_rule_prb_output.txt."
