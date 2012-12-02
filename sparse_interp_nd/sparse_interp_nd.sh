#!/bin/bash
#
mkdir temp
cd temp
rm *
~/bin/$ARCH/f90split ../sparse_interp_nd.f90
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
ar qc libsparse_interp_nd.a *.o
rm *.o
#
mv libsparse_interp_nd.a ~/lib/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/lib/$ARCH/libsparse_interp_nd.a"
