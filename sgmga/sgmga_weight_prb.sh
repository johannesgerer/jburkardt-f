#!/bin/bash
#
gfortran -c -g sgmga_weight_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_weight_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sgmga_weight_prb.o -L$HOME/lib/$ARCH -lsgmga -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_weight_prb.o"
  exit
fi
rm sgmga_weight_prb.o
#
mv a.out sgmga_weight_prb
./sgmga_weight_prb > sgmga_weight_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_weight_prb"
  exit
fi
rm sgmga_weight_prb
#
echo "Program output written to sgmga_weight_prb_output.txt"
