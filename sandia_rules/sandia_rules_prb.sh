#!/bin/bash
#
gfortran -c -g sandia_rules_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_rules_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sandia_rules_prb.o -L$HOME/lib/$ARCH -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_rules_prb.o"
  exit
fi
rm sandia_rules_prb.o
#
mv a.out sandia_rules_prb
./sandia_rules_prb > sandia_rules_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_rules_prb"
  exit
fi
rm sandia_rules_prb
#
echo "Test program output written to sandia_rules_prb_output.txt."
