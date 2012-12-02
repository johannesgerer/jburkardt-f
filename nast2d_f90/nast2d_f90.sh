#!/bin/bash
#
gfortran -c -g nrtype.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nrtype.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g main.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling main.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g boundary.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling boundary.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g extras.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling extras.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g init.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling init.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g surface.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling surface.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g uvp.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uvp.f90."
  exit
fi
rm compiler.txt
#
gfortran -c -g visual.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling visual.f90."
  exit
fi
rm compiler.txt
#
gfortran main.o boundary.o extras.o init.o surface.o uvp.o visual.o
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the nast2d object files."
  exit
fi
rm nrtype.mod
rm *.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/nast2d_f90
#
echo "Executable installed as ~/bin/$ARCH/nast2d_f90."
