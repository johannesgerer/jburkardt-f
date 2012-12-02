#!/bin/bash
#
gfortran -c -g tec_to_obj2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tec_to_obj2.f90"
  exit
fi
rm compiler.txt
#
gfortran tec_to_obj2.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tec_to_obj2.o"
  exit
fi
#
rm tec_to_obj2.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/tec_to_obj2
#
echo "Executable installed as ~/bin/$ARCH/tec_to_obj2"
