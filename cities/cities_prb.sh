#!/bin/bash
#
gfortran -c -g cities_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cities_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cities_prb.o -L$HOME/lib/$ARCH -lcities
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cities_prb.o"
  exit
fi
rm cities_prb.o
#
#  Copy the necessary data files.
#
cp ../../datasets/cities/wg22_main.txt .
cp ../../datasets/cities/wg22_xy.txt .
cp ../../datasets/cities/usca312_main.txt .
cp ../../datasets/cities/usca312_dms.txt .
cp ../../datasets/cities/uscap_ll.txt .
#
#  Run the program.
#
mv a.out cities_prb
./cities_prb > cities_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cities_prb"
  exit
fi
rm cities_prb
#
#  Discard the data files.
#
rm wg22_main.txt
rm wg22_xy.txt
rm usca312_dms.txt
rm usca312_dist.txt
rm usca312_main.txt
rm usca312_xy.txt
rm uscap_ll.txt
#
echo "Program output written to cities_prb_output.txt"
