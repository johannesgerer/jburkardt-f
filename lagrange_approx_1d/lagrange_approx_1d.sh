#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../lagrange_approx_1d.f90
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
ar qc liblagrange_approx_1d.a *.o
rm *.o
#
mv liblagrange_approx_1d.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/liblagrange_approx_1d.a"
