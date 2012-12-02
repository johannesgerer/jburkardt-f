#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../qw_golub_welsch.f90
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
ar qc libqw_golub_welsch.a *.o
rm *.o
#
mv libqw_golub_welsch.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libqw_golub_welsch.a"
