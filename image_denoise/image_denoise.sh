#!/bin/bash
#
mkdir temp
cd temp
rm *
f90split ../image_denoise.f90
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
ar qc libimage_denoise.a *.o
rm *.o
#
mv libimage_denoise.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libimage_denoise.a"
