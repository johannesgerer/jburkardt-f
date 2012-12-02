#!/bin/bash
#
gfortran -c pce_burgers.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pce_burgers.f90"
  exit
fi
rm compiler.txt
#
gfortran pce_burgers.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pce_burgers.o"
  exit
fi
rm pce_burgers.o
#
mv a.out ~/bin/$ARCH/pce_burgers
#
echo "Executable installed as ~/bin/$ARCH/pce_burgers"
