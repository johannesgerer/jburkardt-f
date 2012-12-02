#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../random_data.f90
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
ar qc librandom_data.a *.o
rm *.o
#
mv librandom_data.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/librandom_data.a"
