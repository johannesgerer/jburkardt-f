#!/bin/bash
#
gfortran -c -g dunavant_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dunavant_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran dunavant_prb.o -L$HOME/lib/$ARCH -ldunavant
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dunavant_prb.o"
  exit
fi
rm dunavant_prb.o
#
mv a.out dunavant_prb
./dunavant_prb > dunavant_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dunavant_prb"
  exit
fi
rm dunavant_prb
#
echo "Test program output written to dunavant_prb_output.txt."
#
if [ -e dunavant_rule_01.eps ]; then
  convert dunavant_rule_01.eps dunavant_rule_01.png
  rm dunavant_rule_01.eps
fi
#
if [ -e dunavant_rule_02.eps ]; then
  convert dunavant_rule_02.eps dunavant_rule_02.png
  rm dunavant_rule_02.eps
fi
#
if [ -e dunavant_rule_03.eps ]; then
  convert dunavant_rule_03.eps dunavant_rule_03.png
  rm dunavant_rule_03.eps
fi
#
if [ -e dunavant_rule_04.eps ]; then
  convert dunavant_rule_04.eps dunavant_rule_04.png
  rm dunavant_rule_04.eps
fi
#
if [ -e dunavant_rule_05.eps ]; then
  convert dunavant_rule_05.eps dunavant_rule_05.png
  rm dunavant_rule_05.eps
fi
#
if [ -e dunavant_rule_06.eps ]; then
  convert dunavant_rule_06.eps dunavant_rule_06.png
  rm dunavant_rule_06.eps
fi
#
if [ -e dunavant_rule_07.eps ]; then
  convert dunavant_rule_07.eps dunavant_rule_07.png
  rm dunavant_rule_07.eps
fi
#
if [ -e dunavant_rule_08.eps ]; then
  convert dunavant_rule_08.eps dunavant_rule_08.png
  rm dunavant_rule_08.eps
fi
#
if [ -e dunavant_rule_09.eps ]; then
  convert dunavant_rule_09.eps dunavant_rule_09.png
  rm dunavant_rule_09.eps
fi
#
if [ -e dunavant_rule_10.eps ]; then
  convert dunavant_rule_10.eps dunavant_rule_10.png
  rm dunavant_rule_10.eps
fi
#
echo "Images converted from EPS to PNG format."
