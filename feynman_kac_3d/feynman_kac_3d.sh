#!/bin/bash
#
gfortran -c -g feynman_kac_3d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling feynman_kac_3d.f90"
  exit
fi
rm compiler.txt
#
gfortran feynman_kac_3d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading feynman_kac_3d.o"
  exit
fi
rm feynman_kac_3d.o
#
chmod ugo+x a.out
mv a.out feynman_kac_3d
./feynman_kac_3d > feynman_kac_3d_output.txt
#
rm feynman_kac_3d
#
echo "Program output written to feynman_kac_3d_output.txt"
