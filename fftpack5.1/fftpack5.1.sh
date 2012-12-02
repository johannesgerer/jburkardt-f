#!/bin/bash
#
mkdir temp
cd temp
f90split ../fftpack5.1.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c -O3 $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f90
#
ar cr libfftpack5.1.a *.o
rm *.o
#
mv libfftpack5.1.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libfftpack5.1.a."
