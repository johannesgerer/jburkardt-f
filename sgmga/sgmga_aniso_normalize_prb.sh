#!/bin/bash
#
gfortran -c -g sgmga_aniso_normalize_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_aniso_normalize_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sgmga_aniso_normalize_prb.o -L$HOME/lib/$ARCH -lsgmga -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_aniso_normalize_prb.o"
  exit
fi
rm sgmga_aniso_normalize_prb.o
#
mv a.out sgmga_aniso_normalize_prb
./sgmga_aniso_normalize_prb > sgmga_aniso_normalize_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_aniso_normalize_prb"
  exit
fi
rm sgmga_aniso_normalize_prb
#
echo "Program output written to sgmga_aniso_normalize_prb_output.txt"

