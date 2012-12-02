#!/bin/bash
#
gfortran -c -g triangulation_q2l.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_q2l.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_q2l.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_q2l.o"
  exit
fi
#
rm triangulation_q2l.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_q2l
#
echo "Executable installed as ~/bin/$ARCH/triangulation_q2l"
