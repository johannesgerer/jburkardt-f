#!/bin/bash
#
gfortran -c -g fekete_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fekete_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran fekete_prb.o -L$HOME/lib/$ARCH -lfekete
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fekete_prb.o"
  exit
fi
rm fekete_prb.o
#
mv a.out fekete_prb
./fekete_prb > fekete_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fekete_prb"
  exit
fi
rm fekete_prb
#
echo "Test program output written to fekete_prb_output.txt."
#
if [ -e fekete_rule_1.eps ]; then
  convert fekete_rule_1.eps fekete_rule_1.png
  rm fekete_rule_1.eps
fi
#
if [ -e fekete_rule_2.eps ]; then
  convert fekete_rule_2.eps fekete_rule_2.png
  rm fekete_rule_2.eps
fi
#
if [ -e fekete_rule_3.eps ]; then
  convert fekete_rule_3.eps fekete_rule_3.png
  rm fekete_rule_3.eps
fi
#
if [ -e fekete_rule_4.eps ]; then
  convert fekete_rule_4.eps fekete_rule_4.png
  rm fekete_rule_4.eps
fi
#
if [ -e fekete_rule_5.eps ]; then
  convert fekete_rule_5.eps fekete_rule_5.png
  rm fekete_rule_5.eps
fi
#
if [ -e fekete_rule_6.eps ]; then
  convert fekete_rule_6.eps fekete_rule_6.png
  rm fekete_rule_6.eps
fi
#
if [ -e fekete_rule_7.eps ]; then
  convert fekete_rule_7.eps fekete_rule_7.png
  rm fekete_rule_7.eps
fi
#
echo "EPS output converted to PNG format."
