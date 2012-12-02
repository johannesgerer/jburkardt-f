#!/bin/bash
#
gfortran -c -g sgmga_size_table.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_size_table.f90"
  exit
fi
rm compiler.txt
#
gfortran sgmga_size_table.o -L$HOME/lib/$ARCH -lsgmga -lsandia_rules
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_size_table.o"
  exit
fi
rm sgmga_size_table.o
#
mv a.out sgmga_size_table
./sgmga_size_table > sgmga_size_table_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_size_table"
  exit
fi
rm sgmga_size_table
#
echo "Program output written to sgmga_size_table_output.txt"
