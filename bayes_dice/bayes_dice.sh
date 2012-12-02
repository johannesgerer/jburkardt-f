#!/bin/bash
#
gfortran -c bayes_dice.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bayes_dice.f90"
  exit
fi
rm compiler.txt
#
gfortran bayes_dice.o
if [ $? -ne 0 ]; then
  echo "Errors loading bayes_dice.o"
  exit
fi
rm bayes_dice.o
#
mv a.out ~/bin/$ARCH/bayes_dice
#
echo "Executable installed as ~/bin/$ARCH/bayes_dice"
