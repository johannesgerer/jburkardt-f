#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../partial_digest.f90
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
ar qc libpartial_digest.a *.o
rm *.o
#
mv libpartial_digest.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libpartial_digest.a"
