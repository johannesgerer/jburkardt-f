#!/bin/bash
#
gfortran -c -g sgmga_unique_index_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_unique_index_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sgmga_unique_index_prb.o -L$HOME/lib/$ARCH -lsgmga -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_unique_index_prb.o"
  exit
fi
rm sgmga_unique_index_prb.o
#
mv a.out sgmga_unique_index_prb
./sgmga_unique_index_prb > sgmga_unique_index_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_unique_index_prb"
  exit
fi
rm sgmga_unique_index_prb
#
echo "Program output written to sgmga_unique_index_prb_output.txt"
