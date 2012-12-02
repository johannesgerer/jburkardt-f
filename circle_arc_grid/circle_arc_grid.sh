#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../circle_arc_grid.f90
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
ar qc libcircle_arc_grid.a *.o
rm *.o
#
mv libcircle_arc_grid.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libcircle_arc_grid.a"
