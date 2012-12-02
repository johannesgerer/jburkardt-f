#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../toms291.f90
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
ar qc libtoms291.a *.o
rm *.o
#
mv libtoms291.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libtoms291.a"
