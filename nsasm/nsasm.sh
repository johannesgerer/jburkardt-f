#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../nsasm.f90
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
ar qc libnsasm.a *.o
rm *.o
#
mv libnsasm.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libnsasm.a"
