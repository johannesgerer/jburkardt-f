#!/bin/bash
#
mkdir temp
cd temp
rm *
cp ../mpi_stubs_f90.h .
~/bin/$ARCH/f90split ../mpi_stubs.f90
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
ar qc libmpi_stubs.a *.o
rm *.o
#
mv mpi_stubs_f90.h ~/include
mv libmpi_stubs.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libmpi_stubs.a"
echo "Include file installed as ~/include/mpi_stubs_f90.h"

