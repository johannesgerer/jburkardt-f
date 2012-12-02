#!/bin/bash
#
gfortran -c -g table_record_match.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling table_record_match.f90"
  exit
fi
rm compiler.txt
#
gfortran table_record_match.o -L$HOME/lib/$ARCH
if [ $? -ne 0 ]; then
  echo "Error loading table_record_match.o"
  exit
fi
rm table_record_match.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/table_record_match
#
echo "Executable installed as ~/bin/$ARCH/table_record_match"
