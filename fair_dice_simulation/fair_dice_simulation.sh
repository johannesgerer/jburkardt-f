#!/bin/bash
#
gfortran -c -g fair_dice_simulation.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fair_dice_simulation.f90"
  exit
fi
rm compiler.txt
#
gfortran fair_dice_simulation.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fair_dice_simulation.o"
  exit
fi
rm fair_dice_simulation.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/fair_dice_simulation
#
echo "Executable installed as ~/bin/$ARCH/fair_dice_simulation"
