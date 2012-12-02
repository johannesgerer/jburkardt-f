#!/bin/bash
#
gfortran -c -g chrpak_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chrpak_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran chrpak_prb.o -L$HOME/lib/$ARCH -lchrpak
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chrpak_prb.o"
  exit
fi
rm chrpak_prb.o
#
mv a.out chrpak_prb
./chrpak_prb > chrpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chrpak_prb"
  exit
fi
rm chrpak_prb
#
echo "Test program output written to chrpak_prb_output.txt."
