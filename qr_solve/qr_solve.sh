#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../qr_solve.f90
#
for FILE in `ls -1 *.f90`
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
ar qc libqr_solve.a *.o
rm *.o
#
mv libqr_solve.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libqr_solve.a"
