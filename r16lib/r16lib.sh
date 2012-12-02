#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../r16lib.f90
#
for FILE in `ls -1 *.f90`;
do
  g95 -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f90
#
ar qc libr16lib.a *.o
rm *.o
#
mv libr16lib.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libr16lib.a"
