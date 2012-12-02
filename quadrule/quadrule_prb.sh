#!/bin/bash
#
gfortran -c -g quadrule_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrule_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran quadrule_prb.o -L$HOME/lib/$ARCH -lquadrule
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrule_prb.o"
  exit
fi
rm quadrule_prb.o
#
mv a.out quadrule_prb
./quadrule_prb > quadrule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadrule_prb"
  exit
fi
rm quadrule_prb
#
echo "Test program output written to quadrule_prb_output.txt."
