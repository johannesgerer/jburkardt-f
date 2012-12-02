#!/bin/bash
#
gfortran -c -g sgmga_size_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_size_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sgmga_size_prb.o -L$HOME/lib/$ARCH -lsgmga -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_size_prb.o"
  exit
fi
rm sgmga_size_prb.o
#
mv a.out sgmga_size_prb
./sgmga_size_prb > sgmga_size_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_size_prb"
  exit
fi
rm sgmga_size_prb
#
echo "Program output written to sgmga_size_prb_output.txt"
