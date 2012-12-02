#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../van_der_corput.f90
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
ar qc libvan_der_corput.a *.o
rm *.o
#
mv libvan_der_corput.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libvan_der_corput.a"
