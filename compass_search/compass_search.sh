#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../compass_search.f90
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
ar qc libcompass_search.a *.o
rm *.o
#
mv libcompass_search.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libcompass_search.a"
