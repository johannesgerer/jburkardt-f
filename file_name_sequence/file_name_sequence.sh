#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../file_name_sequence.f90
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
ar qc libfile_name_sequence.a *.o
rm *.o
#
mv libfile_name_sequence.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libfile_name_sequence.a"
