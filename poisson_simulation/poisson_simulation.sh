#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../poisson_simulation.f90
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
ar qc libpoisson_simulation.a *.o
rm *.o
#
mv libpoisson_simulation.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libpoisson_simulation.a"
