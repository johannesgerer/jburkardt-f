#!/bin/bash
#
cp omp_lib.h /$HOME/include
cp omp_lib_kinds.h /$HOME/include
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../openmp_stubs.f90
cp ../omp_lib.h .
cp ../omp_lib_kinds.h .
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f90
#
ar qc libopenmp_stubs.a *.o
rm *.o
rm *.h
#
mv libopenmp_stubs.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libopenmp_stubs.a"
