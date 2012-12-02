#!/bin/bash
#
gfortran -c -g apportionment_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling apportionment_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran apportionment_prb.o -L$HOME/lib/$ARCH -lapportionment
if [ $? -ne 0 ]; then
  echo "Errors linking and loading apportionment_prb.o"
  exit
fi
rm apportionment_prb.o
#
mv a.out apportionment_prb
./apportionment_prb > apportionment_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running apportionment_prb"
  exit
fi
rm apportionment_prb
#
echo "Test program output written to apportionment_prb_output.txt."
