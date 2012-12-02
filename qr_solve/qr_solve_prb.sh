#!/bin/bash
#
gfortran -c -g qr_solve_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qr_solve_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran qr_solve_prb.o -L$HOME/lib/$ARCH -lqr_solve -ltest_ls -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qr_solve_prb.o"
  exit
fi
rm qr_solve_prb.o
#
mv a.out qr_solve_prb
./qr_solve_prb > qr_solve_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qr_solve_prb"
  exit
fi
rm qr_solve_prb
#
echo "Program output written to qr_solve_prb_output.txt"
