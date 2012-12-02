#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../ice_io.f90
#
for FILE in `ls -1 *.f90`;
do
  g95 -c -g -I /usr/local/include $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f90
#
ar qc libice_io.a *.o
rm *.o
#
mv libice_io.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libice_io.a"
