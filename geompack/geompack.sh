#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../geompack.f90
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
ar qc libgeompack.a *.o
rm *.o
#
mv libgeompack.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libgeompack.a"
