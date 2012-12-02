#!/bin/bash
#
gfortran -c -g pyramid_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran pyramid_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_rule.o"
  exit
fi
rm pyramid_rule.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/pyramid_rule
#
echo "Executable installed as ~/bin/$ARCH/pyramid_rule"
