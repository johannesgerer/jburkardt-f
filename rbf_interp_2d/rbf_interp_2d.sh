#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../rbf_interp_2d.f90
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
ar qc librbf_interp_2d.a *.o
rm *.o
#
mv librbf_interp_2d.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/librbf_interp_2d.a"
