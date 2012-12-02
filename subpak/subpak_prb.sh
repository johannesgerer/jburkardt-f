#!/bin/bash
#
gfortran -c -g subpak_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subpak_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran subpak_prb.o -L$HOME/lib/$ARCH -lsubpak
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subpak_prb.o"
  exit
fi
rm subpak_prb.o
#
mv a.out subpak_prb
./subpak_prb > subpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subpak_prb"
  exit
fi
rm subpak_prb
#
echo "Test program output written to subpak_prb_output.txt."
