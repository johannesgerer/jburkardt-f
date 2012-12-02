#!/bin/bash
#
gfortran -c -g sgmga_point_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_point_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sgmga_point_prb.o -L$HOME/lib/$ARCH -lsgmga -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_point_prb.o"
  exit
fi
rm sgmga_point_prb.o
#
mv a.out sgmga_point_prb
./sgmga_point_prb > sgmga_point_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_point_prb"
  exit
fi
rm sgmga_point_prb
#
echo "Program output written to sgmga_point_prb_output.txt"
