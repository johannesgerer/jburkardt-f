#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../linpack_d.f90
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
ar qc liblinpack_d.a *.o
rm *.o
#
mv liblinpack_d.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/liblinpack_d.a"
