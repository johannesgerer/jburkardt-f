#!/bin/bash
#
gfortran -c -g chebyshev_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran chebyshev_prb.o -L$HOME/lib/$ARCH -lchebyshev
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_prb.o"
  exit
fi
rm chebyshev_prb.o
#
mv a.out chebyshev_prb
./chebyshev_prb > chebyshev_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_prb"
  exit
fi
rm chebyshev_prb
#
echo "Program output written to chebyshev_prb_output.txt"
