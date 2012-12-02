#!/bin/bash
#
gfortran -c -g latin_random_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_random_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran latin_random_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_random_dataset.o + liblatin_random.a"
  exit
fi
rm latin_random_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/latin_random_dataset
#
echo "Executable installed as ~/bin/$ARCH/latin_random_dataset"
