#!/bin/bash
#
gfortran -c -g chebyshev_interp_1d_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_interp_1d_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran chebyshev_interp_1d_prb.o -L$HOME/lib/$ARCH -lchebyshev_interp_1d \
  -ltest_interp -lqr_solve -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_interp_1d_prb.o"
  exit
fi
rm chebyshev_interp_1d_prb.o
#
mv a.out chebyshev_interp_1d_prb
./chebyshev_interp_1d_prb > chebyshev_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_interp_1d_prb"
  exit
fi
rm chebyshev_interp_1d_prb
#
echo "Test program output written to chebyshev_interp_1d_prb_output.txt."
