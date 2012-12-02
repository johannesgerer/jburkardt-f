#!/bin/bash
#
gfortran -c -g movie_data_reformat.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling movie_data_reformat.f90"
  exit
fi
rm compiler.txt
#
gfortran movie_data_reformat.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading movie_data_reformat.o"
  exit
fi
rm movie_data_reformat.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/movie_data_reformat
#
echo "Executable installed as ~/bin/$ARCH/movie_data_reformat"
