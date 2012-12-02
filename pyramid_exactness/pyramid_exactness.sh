#!/bin/bash
#
gfortran -c -g pyramid_exactness.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_exactness.f90"
  exit
fi
rm compiler.txt
#
gfortran pyramid_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_exactness.o"
  exit
fi
rm pyramid_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/pyramid_exactness
#
echo "Executable installed as ~/bin/$ARCH/pyramid_exactness"
