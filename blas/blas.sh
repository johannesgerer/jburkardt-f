#!/bin/bash
#
mkdir temp
cd temp
f90split ../../blas1_c/blas1_c.f90
f90split ../../blas1_d/blas1_d.f90
f90split ../../blas1_s/blas1_s.f90
f90split ../../blas1_z/blas1_z.f90
f90split ../../blas2/blas2.f90
f90split ../../blas3/blas3.f90
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
#  Create the library.
#
ar qc libblas.a *.o
rm *.o
#
mv libblas.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libblas.a."
