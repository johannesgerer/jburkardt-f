#!/bin/bash
#
gfortran -c -g set_theory_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling set_theory_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran set_theory_prb.o -L$HOME/lib/$ARCH -lset_theory
if [ $? -ne 0 ]; then
  echo "Errors linking and loading set_theory_prb.o"
  exit
fi
rm set_theory_prb.o
#
mv a.out set_theory_prb
./set_theory_prb > set_theory_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running set_theory_prb"
  exit
fi
rm set_theory_prb
#
echo "Test program output written to set_theory_prb_output.txt."
