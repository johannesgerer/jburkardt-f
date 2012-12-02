#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../arpack.f90
~/bin/$ARCH/f90split ../blas.f90
~/bin/$ARCH/f90split ../lapack.f90
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
ar qc libarpack.a *.o
rm *.o
#
mv libarpack.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libarpack.a"
