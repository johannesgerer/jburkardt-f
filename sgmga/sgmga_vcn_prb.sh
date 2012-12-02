#!/bin/bash
#
gfortran -c -g sgmga_vcn_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_vcn_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sgmga_vcn_prb.o -L$HOME/lib/$ARCH -lsgmga -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_vcn_prb.o"
  exit
fi
rm sgmga_vcn_prb.o
#
mv a.out sgmga_vcn_prb
./sgmga_vcn_prb > sgmga_vcn_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_vcn_prb"
  exit
fi
rm sgmga_vcn_prb
#
echo "Program output written to sgmga_vcn_prb_output.txt"

