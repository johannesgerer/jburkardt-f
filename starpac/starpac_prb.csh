#!/bin/csh
#
F90 -c -g starpac_prb.f90 >& compiler.txt
if ( $status != 0 ) then
  echo "Errors compiling starpac_prb.f"
  exit
endif
rm compiler.txt
#
F90 starpac_prb.o -L$HOME/lib/$ARCH -lstarpac
if ( $status != 0 ) then
  echo "Errors linking and loading starpac_prb.o"
  exit
endif
rm starpac_prb.o
#
mv a.out starpac_prb
./starpac_prb > starpac_prb_output.txt
if ( $status != 0 ) then
  echo "Errors running starpac_prb"
  exit
endif
rm starpac_prb
#
echo "Program output written to starpac_prb_output.txt"
