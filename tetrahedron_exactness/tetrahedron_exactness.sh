#!/bin/bash
#
gfortran -c -g tetrahedron_exactness.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_exactness.f90"
  exit
fi
rm compiler.txt
#
gfortran tetrahedron_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_exactness.o"
  exit
fi
rm tetrahedron_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/tetrahedron_exactness
#
echo "Program installed as ~/bin/$ARCH/tetrahedron_exactness"
