#!/bin/bash
#
gfortran -c -g gm_rule_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gm_rule_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran gm_rule_prb.o -L$HOME/lib/$ARCH -lgm_rule
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gm_rule_prb.o"
  exit
fi
rm gm_rule_prb.o
#
mv a.out gm_rule_prb
./gm_rule_prb > gm_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gm_rule_prb"
  exit
fi
rm gm_rule_prb
#
echo "Test program output written to gm_rule_prb_output.txt."
