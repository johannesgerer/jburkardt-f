#!/bin/bash
#
gfortran -c -g u_to_f_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling u_to_f_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran u_to_f_prb.o -L$HOME/lib/$ARCH -lrejoin
if [ $? -ne 0 ]; then
  echo "Errors linking and loading u_to_f_prb.o"
  exit
fi
rm u_to_f_prb.o
#
mv a.out u_to_f_prb
./u_to_f_prb > u_to_f_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running u_to_f_prb"
  exit
fi
rm u_to_f_prb
#
echo "Test program output written to u_to_f_prb_output.txt."
