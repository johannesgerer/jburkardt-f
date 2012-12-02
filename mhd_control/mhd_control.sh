#!/bin/bash
#
gfortran -c -g mhd_control.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mhd_control.f90"
  exit
fi
rm compiler.txt
#
gfortran mhd_control.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mhd_control.o"
  exit
fi
rm mhd_control.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/mhd_control
#
echo "A new version of mhd_control has been created."
