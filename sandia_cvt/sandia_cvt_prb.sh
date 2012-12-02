#!/bin/bash
#
gfortran -c -g sandia_cvt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_cvt_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sandia_cvt_prb.o -L$HOME/lib/$ARCH -lsandia_cvt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_cvt_prb.o"
  exit
fi
rm sandia_cvt_prb.o
#
mv a.out sandia_cvt_prb
./sandia_cvt_prb > sandia_cvt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_cvt_prb"
  exit
fi
rm sandia_cvt_prb
#
echo "Test program output written to sandia_cvt_prb_output.txt."
