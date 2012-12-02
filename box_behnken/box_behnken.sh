#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../box_behnken.f90
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
ar qc libbox_behnken.a *.o
rm *.o
#
mv libbox_behnken.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libbox_behnken.a"
