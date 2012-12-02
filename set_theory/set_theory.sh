#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../set_theory.f90
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
ar qc libset_theory.a *.o
rm *.o
#
mv libset_theory.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libset_theory.a"
