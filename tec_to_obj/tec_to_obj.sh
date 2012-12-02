#!/bin/bash
#
gfortran -c -g tec_to_obj.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tec_to_obj.f90"
  exit
fi
rm compiler.txt
#
gfortran tec_to_obj.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tec_to_obj.o"
  exit
fi
#
rm tec_to_obj.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/tec_to_obj
#
echo "Executable installed as ~/bin/$ARCH/tec_to_obj"
