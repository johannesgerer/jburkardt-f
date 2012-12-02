#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../sparse_grid_hw.f90
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
ar qc libsparse_grid_hw.a *.o
rm *.o
#
mv libsparse_grid_hw.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libsparse_grid_hw.a"
